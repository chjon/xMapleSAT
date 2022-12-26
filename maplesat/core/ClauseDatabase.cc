#include "core/ClauseDatabase.h"
#include "core/Solver.h"
#include "mtl/Sort.h"

using namespace Minisat;

static const char* _cat = "CORE";
static DoubleOption opt_garbage_frac (_cat, "gc-frac", "The fraction of wasted memory allowed before a garbage collection is triggered", 0.20, DoubleRange(0, false, HUGE_VAL, false));

ClauseDatabase::ClauseDatabase(Solver* s)
    : remove_satisfied(true)
    , garbage_frac (opt_garbage_frac)

    // Statistics
    , clauses_literals(0)
    , learnts_literals(0)

    // Solver references
    , variableDatabase(s->variableDatabase)
    , ca(s->ca)
    , assignmentTrail(s->assignmentTrail)
    , unitPropagator(s->unitPropagator)
    , branchingHeuristicManager(s->branchingHeuristicManager)
    , solver(s)
{}

void ClauseDatabase::removeClause(CRef cr) {
    Clause& c = ca[cr];
    unitPropagator.detachClause(c, cr);

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
    solver->relocAll(to);

    if (solver->verbosity >= 2)
        printf("|  Garbage collection:   %12d bytes => %12d bytes             |\n", 
            ca.size() * ClauseAllocator::Unit_Size, to.size() * ClauseAllocator::Unit_Size);
    
    // Transfer ownership of memory
    to.moveTo(ca);
}

/*_________________________________________________________________________________________________
|
|  reduceDB : ()  ->  [void]
|  
|  Description:
|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
|________________________________________________________________________________________________@*/
struct reduceDB_lt { 
    ClauseAllocator& ca;
#if LBD_BASED_CLAUSE_DELETION
    const vec<double>& activity;
    reduceDB_lt(ClauseAllocator& ca_,const vec<double>& activity_) : ca(ca_), activity(activity_) {}
#else
    reduceDB_lt(ClauseAllocator& ca_) : ca(ca_) {}
#endif
    bool operator () (CRef x, CRef y) { 
#if LBD_BASED_CLAUSE_DELETION
        return ca[x].activity() > ca[y].activity();
#else
        return ca[x].size() > 2 && (ca[y].size() == 2 || ca[x].activity() < ca[y].activity());
#endif
    }
};

void ClauseDatabase::reduceDB() {
    // Sort clauses by activity
#if LBD_BASED_CLAUSE_DELETION
    sort(learnts, reduceDB_lt(ca, branchingHeuristicManager.getActivityVSIDS()));
#else
    double extra_lim = cla_inc / learnts.size(); // Remove any clause below this activity
    sort(learnts, reduceDB_lt(ca));
#endif

    // Don't delete binary or locked clauses. From the rest, delete clauses from the first half
    // and clauses with activity smaller than 'extra_lim':
    int i, j;
    for (i = j = 0; i < learnts.size(); i++){
        Clause& c = ca[learnts[i]];
#if LBD_BASED_CLAUSE_DELETION
        if (c.activity() > 2 && !assignmentTrail.locked(c) && i < learnts.size() / 2)
#else
        if (c.size() > 2 && !assignmentTrail.locked(c) && (i < learnts.size() / 2 || c.activity() < extra_lim))
#endif
            removeClause(learnts[i]);
        else
            learnts[j++] = learnts[i];
    }
    learnts.shrink(i - j);

    // Perform garbage collection if needed
    checkGarbage(garbage_frac);
}

void ClauseDatabase::removeSatisfied(vec<CRef>& cs) {
    int i, j;
    for (i = j = 0; i < cs.size(); i++) {
        Clause& c = ca[cs[i]];
        if (variableDatabase.satisfied(c))
            removeClause(cs[i]);
        else
            cs[j++] = cs[i];
    }
    cs.shrink(i - j);
}

//=================================================================================================
// Writing CNF to DIMACS:
// 
// FIXME: this needs to be rewritten completely.

static Var mapVar(Var x, vec<Var>& map, Var& max) {
    if (map.size() <= x || map[x] == -1){
        map.growTo(x+1, -1);
        map[x] = max++;
    }
    return map[x];
}


void ClauseDatabase::toDimacs(FILE* f, Clause& c, vec<Var>& map, Var& max) {
    if (variableDatabase.satisfied(c)) return;

    for (int i = 0; i < c.size(); i++)
        if (variableDatabase.value(c[i]) != l_False)
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
    if (!solver->ok){
        fprintf(f, "p cnf 1 2\n1 0\n-1 0\n");
        return; }

    vec<Var> map; Var max = 0;

    // Cannot use removeClauses here because it is not safe
    // to deallocate them at this point. Could be improved.
    int cnt = 0;
    for (int i = 0; i < clauses.size(); i++)
        if (!variableDatabase.satisfied(ca[clauses[i]]))
            cnt++;
        
    for (int i = 0; i < clauses.size(); i++)
        if (!variableDatabase.satisfied(ca[clauses[i]])){
            Clause& c = ca[clauses[i]];
            for (int j = 0; j < c.size(); j++)
                if (variableDatabase.value(c[j]) != l_False)
                    mapVar(var(c[j]), map, max);
        }

    // Assumptions are added as unit clauses:
    cnt += solver->assumptions.size();

    fprintf(f, "p cnf %d %d\n", max, cnt);

    for (int i = 0; i < solver->assumptions.size(); i++){
        assert(variableDatabase.value(solver->assumptions[i]) != l_False);
        fprintf(f, "%s%d 0\n", sign(solver->assumptions[i]) ? "-" : "", mapVar(var(solver->assumptions[i]), map, max)+1);
    }

    for (int i = 0; i < clauses.size(); i++)
        toDimacs(f, ca[clauses[i]], map, max);

    if (solver->verbosity > 0)
        printf("Wrote %d clauses with %d variables.\n", cnt, max);
}
