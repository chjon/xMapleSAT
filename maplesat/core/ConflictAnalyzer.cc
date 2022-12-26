#include "core/ConflictAnalyzer.h"
#include "core/Solver.h"

using namespace Minisat;

static const char* _cat = "CORE";

static IntOption opt_ccmin_mode (_cat, "ccmin-mode", "Controls conflict clause minimization (0=none, 1=basic, 2=deep)", 2, IntRange(0, 2));

ConflictAnalyzer::ConflictAnalyzer(Solver* s)
    /////////////
    // Parameters

    : ccmin_mode(opt_ccmin_mode)

    /////////////
    // Statistics

    , max_literals(0)
    , tot_literals(0)

    ////////////////////
    // Solver references

    , assignmentTrail(s->assignmentTrail)
    , branchingHeuristicManager(s->branchingHeuristicManager)
    , ca(s->ca)
    , solver(s)
{}

// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool ConflictAnalyzer::litRedundant(Lit p, uint32_t abstract_levels) {
    analyze_stack.clear(); analyze_stack.push(p);
    int top = solver->analyze_toclear.size();
    while (analyze_stack.size() > 0){
        assert(assignmentTrail.reason(var(analyze_stack.last())) != CRef_Undef);
        Clause& c = ca[assignmentTrail.reason(var(analyze_stack.last()))]; analyze_stack.pop();

        for (int i = 1; i < c.size(); i++){
            Lit p = c[i];
            if (solver->seen[var(p)] || assignmentTrail.level(var(p)) == 0) continue;

            if (assignmentTrail.reason(var(p)) != CRef_Undef && (assignmentTrail.abstractLevel(var(p)) & abstract_levels) != 0){
                solver->seen[var(p)] = 1;
                analyze_stack.push(p);
                solver->analyze_toclear.push(p);
            } else {
                for (int j = top; j < solver->analyze_toclear.size(); j++)
                    solver->seen[var(solver->analyze_toclear[j])] = 0;
                solver->analyze_toclear.shrink(solver->analyze_toclear.size() - top);
                return false;
            }
        }
    }

    return true;
}

inline void ConflictAnalyzer::getFirstUIPClause(CRef confl, vec<Lit>& learntClause) {
    Lit p = lit_Undef;

    // Generate conflict clause by backtracking until the first UIP:
    //
    learntClause.push(); // (leave room for the asserting literal)
    int index = assignmentTrail.nAssigns() - 1;

    // Iterate through every clause that participates in the conflict graph
    int pathC = 0;
    do {
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

#if LBD_BASED_CLAUSE_DELETION
        if (c.learnt() && c.activity() > 2)
            c.activity() = solver->lbd(c);
#else
        if (c.learnt())
            claBumpActivity(c);
#endif

        // Iterate through every literal that participates in the conflict graph
        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];
            if (solver->seen[var(q)] || assignmentTrail.level(var(q)) == 0) continue;

            // Mark variable as seen
            branchingHeuristicManager.handleEventLitInConflictGraph(q, solver->conflicts);
            solver->seen[var(q)] = 1;

            // Increment the number of paths if the variable is assigned after the decision level that resulted in conflict
            if (assignmentTrail.level(var(q)) == assignmentTrail.decisionLevel())
                pathC++;
            // Add literals from earlier decision levels to the conflict clause
            else
                learntClause.push(q);
        }
        
        // Select next clause to look at:
        while (!solver->seen[var(assignmentTrail[index--])]);
        p     = assignmentTrail[index+1];
        confl = assignmentTrail.reason(var(p));
        solver->seen[var(p)] = 0;
        pathC--;

    } while (pathC > 0);

    // Add first UIP literal at index 0
    learntClause[0] = ~p;
}

inline void ConflictAnalyzer::simplifyClause(vec<Lit>& learntClause) {
    int i, j;
    if (ccmin_mode == 2) {
        uint32_t abstract_level = 0;
        for (i = 1; i < learntClause.size(); i++)
            abstract_level |= assignmentTrail.abstractLevel(var(learntClause[i])); // (maintain an abstraction of levels involved in conflict)

        for (i = j = 1; i < learntClause.size(); i++)
            if (assignmentTrail.reason(var(learntClause[i])) == CRef_Undef || !litRedundant(learntClause[i], abstract_level))
                learntClause[j++] = learntClause[i];
        
    } else if (ccmin_mode == 1) {
        for (i = j = 1; i < learntClause.size(); i++){
            Var x = var(learntClause[i]);

            if (assignmentTrail.reason(x) == CRef_Undef) {
                learntClause[j++] = learntClause[i];
            } else {
                Clause& c = ca[assignmentTrail.reason(var(learntClause[i]))];
                for (int k = 1; k < c.size(); k++)
                    if (!solver->seen[var(c[k])] && assignmentTrail.level(var(c[k])) > 0){
                        learntClause[j++] = learntClause[i];
                        break; }
            }
        }
    } else {
        i = j = learntClause.size();
    }

    learntClause.shrink(i - j);
}

inline void ConflictAnalyzer::enforceWatcherInvariant(vec<Lit>& learntClause) {
    // Nothing to do for unit clauses
    if (learntClause.size() == 1) return;

    // Find the first literal assigned at the next-highest level:
    int max_i = 1;
    for (int i = 2; i < learntClause.size(); i++) {
        if (assignmentTrail.level(var(learntClause[i])) > assignmentTrail.level(var(learntClause[max_i])))
            max_i = i;
    }

    // Swap-in this literal at index 1:
    std::swap(learntClause[1], learntClause[max_i]);
}

/*_________________________________________________________________________________________________
|
|  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
|  
|  Description:
|    Analyze conflict and produce a reason clause.
|  
|    Pre-conditions:
|      * 'out_learnt' is assumed to be cleared.
|      * Current decision level must be greater than root level.
|  
|    Post-conditions:
|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
|      * If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of the 
|        rest of literals. There may be others from the same level though.
|  
|________________________________________________________________________________________________@*/
void ConflictAnalyzer::analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel) {
    // Generate conflict clause:
    getFirstUIPClause(confl, out_learnt);
    max_literals += out_learnt.size();

    // Simplify conflict clause:
    out_learnt.copyTo(solver->analyze_toclear);
    simplifyClause(out_learnt);
    tot_literals += out_learnt.size();

    // Enforce watcher invariant
    enforceWatcherInvariant(out_learnt);

    // Find backtrack level:
    out_btlevel = (out_learnt.size() == 1) ? 0 : assignmentTrail.level(var(out_learnt[1]));

    // Update data structures for branching heuristics
    branchingHeuristicManager.handleEventLearnedClause(out_learnt);

    // Clear 'seen[]'
    for (int j = 0; j < solver->analyze_toclear.size(); j++)
        solver->seen[var(solver->analyze_toclear[j])] = 0;
}

/*_________________________________________________________________________________________________
|
|  analyzeFinal : (p : Lit)  ->  [void]
|  
|  Description:
|    Specialized analysis procedure to express the final conflict in terms of assumptions.
|    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
|    stores the result in 'out_conflict'.
|________________________________________________________________________________________________@*/
void ConflictAnalyzer::analyzeFinal(Lit p, vec<Lit>& out_conflict) {
    out_conflict.clear();
    out_conflict.push(p);

    if (assignmentTrail.decisionLevel() == 0)
        return;

    solver->seen[var(p)] = 1;

    for (int i = assignmentTrail.nAssigns() - 1; i >= assignmentTrail.indexOfDecisionLevel(1); i--) {
        Var x = var(assignmentTrail[i]);
        if (solver->seen[x]) {
            if (assignmentTrail.reason(x) == CRef_Undef){
                assert(assignmentTrail.level(x) > 0);
                out_conflict.push(~assignmentTrail[i]);
            } else {
                Clause& c = ca[assignmentTrail.reason(x)];
                for (int j = 1; j < c.size(); j++)
                    if (assignmentTrail.level(var(c[j])) > 0)
                        solver->seen[var(c[j])] = 1;
            }
            solver->seen[x] = 0;
        }
    }

    solver->seen[var(p)] = 0;
}