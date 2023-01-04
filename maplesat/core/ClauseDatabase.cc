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

static const char* _cat = "CORE";
static DoubleOption opt_garbage_frac (_cat, "gc-frac", "The fraction of wasted memory allowed before a garbage collection is triggered", 0.20, DoubleRange(0, false, HUGE_VAL, false));

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

ClauseDatabase::ClauseDatabase(Solver& s)
    // Solver references
    : ca(s.ca)
    , assignmentTrail(s.assignmentTrail)
    , unitPropagator(s.unitPropagator)
    , solver(s)

    // Database growth timer
    , learntSizeLimitGrowthTimer    (100)
    
    // Memory management parameters
    , remove_satisfied(true)
    , garbage_frac (opt_garbage_frac)

    // Database growth parameters
    , learntSizeTimerGrowthFactor   (1.5)
#if !RAPID_DELETION
    , learntSizeLimitFactorInitial  (1./3.)
    , learntSizeLimitGrowthFactor   (1.1)
#endif

    // Statistics
    , clauses_literals(0)
    , learnts_literals(0)
{}

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
    if (!solver.ok){
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
        if (!assignmentTrail.satisfied(ca[clauses[i]])){
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
        printf("Wrote %d clauses with %d variables.\n", cnt, max);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

void ClauseDatabase::removeClause(CRef cr) {
    Clause& c = ca[cr];
    unitPropagator.detachClause(cr);

    // Update stats
    if (c.learnt()) learnts_literals -= c.size();
    else            clauses_literals -= c.size();

    // Don't leave pointers to free'd memory!
    assignmentTrail.handleEventClauseDeleted(c);

    // Mark clause as deleted
    c.mark(1); 
    ca.free(cr);
}

void ClauseDatabase::garbageCollect(void) {
    // Initialize the next region to a size corresponding to the estimated utilization degree. This
    // is not precise but should avoid some unnecessary reallocations for the new region:
    ClauseAllocator to(ca.size() - ca.wasted()); 

    // Reloc all clause references
    unitPropagator .relocAll(to);
    assignmentTrail.relocAll(to);
    relocAll(to);

    if (solver.verbosity >= 2)
        printf("|  Garbage collection:   %12d bytes => %12d bytes             |\n", 
            ca.size() * ClauseAllocator::Unit_Size, to.size() * ClauseAllocator::Unit_Size);
    
    // Transfer ownership of memory
    to.moveTo(ca);
}

void ClauseDatabase::handleEventLearntClause(void) {
    if (--learntSizeLimitGrowthTimerCounter != 0) return;

    // Compute the next time the clause database should grow
    learntSizeLimitGrowthTimer       *= learntSizeTimerGrowthFactor;
    learntSizeLimitGrowthTimerCounter = (int)learntSizeLimitGrowthTimer;

#if ! RAPID_DELETION
    // Update the maximum size of the clause database
    maxNumLearnts *= learntSizeLimitGrowthFactor;
#endif

    if (solver.verbosity >= 1)
        printf(
            "| %9d | %7d %8d %8d | %8d %8d %6.0f | %6.3f %% |\n", 
            (int)solver.conflicts, 
            solver.nFreeVars(),
            nClauses(),
            (int)clauses_literals, 
            (int)maxNumLearnts,
            nLearnts(),
            (double)learnts_literals / nLearnts(),
            assignmentTrail.progressEstimate() * 100
        );
}

struct reduceDB_lbdDeletion_lt {
    const ClauseAllocator& ca;
    const vec<double>& activity;
    reduceDB_lbdDeletion_lt(ClauseAllocator& ca_,const vec<double>& activity_)
        : ca(ca_)
        , activity(activity_)
    {}
    bool operator () (CRef x, CRef y) { 
        return ca[x].activity() > ca[y].activity();
    }
};

struct reduceDB_activityDeletion_lt {
    const ClauseAllocator& ca;
    reduceDB_activityDeletion_lt(ClauseAllocator& ca_)
        : ca(ca_)
    {}
    bool operator () (CRef x, CRef y) { 
        return
            ca[x].size() > 2 && (ca[y].size() == 2 ||
            ca[x].activity() < ca[y].activity());
    }
};

void ClauseDatabase::preprocessReduceDB(void) {
    // Sort clauses by activity
#if LBD_BASED_CLAUSE_DELETION
    sort(learnts, reduceDB_lbdDeletion_lt(ca, solver.branchingHeuristicManager.getActivityVSIDS()));
#else
    extra_lim = cla_inc / learnts.size(); // Remove any clause below this activity
    sort(learnts, reduceDB_activityDeletion_lt(ca));
#endif
}