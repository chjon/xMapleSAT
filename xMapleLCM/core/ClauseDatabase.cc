/*******************************************************************************[ClauseDatabase.cc]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010, Niklas Sorensson

MapleSAT_Refactor, based on MapleSAT -- Copyright (c) 2022, Jonathan Chung, Vijay Ganesh, Sam Buss

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include "core/ClauseDatabase.h"
#include "core/Solver.h"
#include "mtl/Sort.h"

using namespace Minisat;

///////////////////////////////////////////////////////////////////////////////////////////////////
// OPTIONS

static const char* _cat = "CORE";

static DoubleOption opt_garbage_frac (_cat, "gc-frac", "The fraction of wasted memory allowed before a garbage collection is triggered", 0.20, DoubleRange(0, false, HUGE_VAL, false));
static DoubleOption opt_clause_decay (_cat, "cla-decay", "The clause activity decay factor", 0.999, DoubleRange(0, false, 1, false));

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

ClauseDatabase::ClauseDatabase(Solver& s)
    // Solver references
    : ca(s.ca)
    , assignmentTrail(s.assignmentTrail)
    , unitPropagator(s.unitPropagator)
    , solver(s)

    // Memory management parameters
    , remove_satisfied(true)
    , garbage_frac (opt_garbage_frac)

    , next_T2_reduce(10000)
    , next_L_reduce (15000)

    // Clause deletion heuristic parameters
    , cla_inc(1)
    , clause_decay(opt_clause_decay)
    , core_lbd_cut(3)

    // Statistics
    , clauses_literals(0)
    , learnts_literals(0)
    , lbd_calls(0)
{}

void ClauseDatabase::removeClause(CRef cr) {
    Clause& c = ca[cr];
    solver.handleEventClauseDeleted(c, cr);

    // Update stats
    if (c.learnt()) learnts_literals -= c.size();
    else            clauses_literals -= c.size();

    c.mark(1);
    ca.free(cr);
}

void ClauseDatabase::attachClause(CRef cr) {
    const Clause& c = ca[cr];
    unitPropagator.attachClause(cr);
    if (c.learnt()) learnts_literals += c.size();
    else            clauses_literals += c.size();
}

void ClauseDatabase::detachClause(CRef cr, bool strict) {
    const Clause& c = ca[cr];
    unitPropagator.detachClause(cr, strict);

    if (c.learnt()) learnts_literals -= c.size();
    else            clauses_literals -= c.size();
}

struct reduceDB_lt { 
    ClauseAllocator& ca;
    reduceDB_lt(ClauseAllocator& ca_) : ca(ca_) {}
    bool operator() (CRef x, CRef y) const { return ca[x].activity() < ca[y].activity(); }
};

void ClauseDatabase::reduceDB() {
    int     i, j;
    //if (local_learnts_dirty) cleanLearnts(learnts_local, LOCAL);
    //local_learnts_dirty = false;

    sort(learnts_local, reduceDB_lt(ca));

    int limit = learnts_local.size() / 2;
    for (i = j = 0; i < learnts_local.size(); i++){
        Clause& c = ca[learnts_local[i]];
        if (c.mark() == LOCAL)
            if (c.removable() && !assignmentTrail.locked(c) && i < limit)
                removeClause(learnts_local[i]);
            else{
                if (!c.removable()) limit++;
                c.removable(true);
                learnts_local[j++] = learnts_local[i]; }
    }
    learnts_local.shrink(i - j);

    checkGarbage();
}

void ClauseDatabase::reduceDB_Tier2() {
    int i, j;
    for (i = j = 0; i < learnts_tier2.size(); i++) {
        Clause& c = ca[learnts_tier2[i]];
        if (c.mark() == TIER2)
            if (!assignmentTrail.locked(c) && c.touched() + 30000 < solver.conflicts){
                learnts_local.push(learnts_tier2[i]);
                c.mark(LOCAL);
                //c.removable(true);
                c.activity() = 0;
                claBumpActivity(c);
            } else
                learnts_tier2[j++] = learnts_tier2[i];
    }
    learnts_tier2.shrink(i - j);
}


void ClauseDatabase::removeSatisfied(vec<CRef>& cs) {
    int i, j;
    for (i = j = 0; i < cs.size(); i++) {
        Clause& c = ca[cs[i]];
        if (assignmentTrail.satisfied(c))
            removeClause(cs[i]);
        else
            cs[j++] = cs[i];
    }
    cs.shrink(i - j);
}

void ClauseDatabase::safeRemoveSatisfied(vec<CRef>& cs, unsigned valid_mark) {
    int i, j;
    for (i = j = 0; i < cs.size(); i++) {
        Clause& c = ca[cs[i]];
        if (c.mark() == valid_mark)
            if (assignmentTrail.satisfied(c))
                removeClause(cs[i]);
            else
                cs[j++] = cs[i];
    }
    cs.shrink(i - j);
}

void ClauseDatabase::removeSatisfied(void) {
    removeSatisfied(learnts_core); // Should clean core first.
    safeRemoveSatisfied(learnts_tier2, TIER2);
    safeRemoveSatisfied(learnts_local, LOCAL);
    if (remove_satisfied)        // Can be turned off.
        removeSatisfied(clauses);
    checkGarbage();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// OUTPUT

// FIXME: this needs to be rewritten completely.

static Var mapVar(Var x, vec<Var>& map, Var& max) {
    if (map.size() <= x || map[x] == -1){
        map.growTo(x+1, -1);
        map[x] = max++;
    }
    return map[x];
}


void ClauseDatabase::toDimacs(FILE* f, Clause& c, vec<Var>& map, Var& max) {
    if (assignmentTrail.satisfied(c)) return;

    for (int i = 0; i < c.size(); i++)
        if (assignmentTrail.value(c[i]) != l_False)
            fprintf(f, "%s%d ", sign(c[i]) ? "-" : "", mapVar(var(c[i]), map, max)+1);
    fprintf(f, "0\n");
}


void ClauseDatabase::toDimacs(const char *file, const vec<Lit>& assumps) {
    FILE* f = fopen(file, "wr");
    if (f == NULL)
        fprintf(stderr, "could not open file %s\n", file), exit(1);
    toDimacs(f, assumps);
    fclose(f);
}


void ClauseDatabase::toDimacs(FILE* f, const vec<Lit>& assumps) {
    // Handle case when solver is in contradictory state:
    if (!solver.okay()){
        fprintf(f, "p cnf 1 2\n1 0\n-1 0\n");
        return; }

    vec<Var> map; Var max = 0;

    // Cannot use removeClauses here because it is not safe
    // to deallocate them at this point. Could be improved.
    int cnt = 0;
    for (int i = 0; i < clauses.size(); i++)
        if (!assignmentTrail.satisfied(ca[clauses[i]]))
            cnt++;

    for (int i = 0; i < clauses.size(); i++)
        if (!assignmentTrail.satisfied(ca[clauses[i]])) {
            Clause& c = ca[clauses[i]];
            for (int j = 0; j < c.size(); j++)
                if (assignmentTrail.value(c[j]) != l_False)
                    mapVar(var(c[j]), map, max);
        }

    // Assumptions are added as unit clauses:
    cnt += solver.assumptions.size();

    fprintf(f, "p cnf %d %d\n", max, cnt);

    for (int i = 0; i < solver.assumptions.size(); i++){
        assert(assignmentTrail.value(solver.assumptions[i]) != l_False);
        fprintf(f, "%s%d 0\n", sign(solver.assumptions[i]) ? "-" : "", mapVar(var(solver.assumptions[i]), map, max)+1);
    }

    for (int i = 0; i < clauses.size(); i++)
        toDimacs(f, ca[clauses[i]], map, max);

    if (solver.verbosity > 0)
        printf("c Wrote %d clauses with %d variables.\n", cnt, max);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

void ClauseDatabase::garbageCollect(void) {
    // Initialize the next region to a size corresponding to the estimated utilization degree. This
    // is not precise but should avoid some unnecessary reallocations for the new region:
    ClauseAllocator to(ca.size() - ca.wasted()); 

    // Reloc all clause references
    solver.relocAll(to);

    if (solver.verbosity >= 2)
        printf("|  Garbage collection:   %12d bytes => %12d bytes             |\n", 
            ca.size() * ClauseAllocator::Unit_Size, to.size() * ClauseAllocator::Unit_Size);
    
    // Transfer ownership of memory
    to.moveTo(ca);
}
