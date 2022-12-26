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
    int pathC = 0;
    Lit p     = lit_Undef;

    // Generate conflict clause:
    //
    out_learnt.push();      // (leave room for the asserting literal)
    int index = assignmentTrail.nAssigns() - 1;

    do{
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

#if LBD_BASED_CLAUSE_DELETION
        if (c.learnt() && c.activity() > 2)
            c.activity() = solver->lbd(c);
#else
        if (c.learnt())
            claBumpActivity(c);
#endif

        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];

            if (!solver->seen[var(q)] && assignmentTrail.level(var(q)) > 0){
                branchingHeuristicManager.handleEventLitInConflictGraph(q, solver->conflicts);
                solver->seen[var(q)] = 1;
                if (assignmentTrail.level(var(q)) >= assignmentTrail.decisionLevel())
                    pathC++;
                else
                    out_learnt.push(q);
            }
        }
        
        // Select next clause to look at:
        while (!solver->seen[var(assignmentTrail[index--])]);
        p     = assignmentTrail[index+1];
        confl = assignmentTrail.reason(var(p));
        solver->seen[var(p)] = 0;
        pathC--;

    } while (pathC > 0);
    out_learnt[0] = ~p;

    // Simplify conflict clause:
    //
    int i, j;
    out_learnt.copyTo(solver->analyze_toclear);
    if (ccmin_mode == 2){
        uint32_t abstract_level = 0;
        for (i = 1; i < out_learnt.size(); i++)
            abstract_level |= assignmentTrail.abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)

        for (i = j = 1; i < out_learnt.size(); i++)
            if (assignmentTrail.reason(var(out_learnt[i])) == CRef_Undef || !litRedundant(out_learnt[i], abstract_level))
                out_learnt[j++] = out_learnt[i];
        
    }else if (ccmin_mode == 1){
        for (i = j = 1; i < out_learnt.size(); i++){
            Var x = var(out_learnt[i]);

            if (assignmentTrail.reason(x) == CRef_Undef)
                out_learnt[j++] = out_learnt[i];
            else{
                Clause& c = ca[assignmentTrail.reason(var(out_learnt[i]))];
                for (int k = 1; k < c.size(); k++)
                    if (!solver->seen[var(c[k])] && assignmentTrail.level(var(c[k])) > 0){
                        out_learnt[j++] = out_learnt[i];
                        break; }
            }
        }
    }else
        i = j = out_learnt.size();

    max_literals += out_learnt.size();
    out_learnt.shrink(i - j);
    tot_literals += out_learnt.size();

    // Find correct backtrack level:
    //
    if (out_learnt.size() == 1)
        out_btlevel = 0;
    else{
        int max_i = 1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 2; i < out_learnt.size(); i++)
            if (assignmentTrail.level(var(out_learnt[i])) > assignmentTrail.level(var(out_learnt[max_i])))
                max_i = i;
        // Swap-in this literal at index 1:
        Lit p             = out_learnt[max_i];
        out_learnt[max_i] = out_learnt[1];
        out_learnt[1]     = p;
        out_btlevel       = assignmentTrail.level(var(p));
    }

    branchingHeuristicManager.handleEventLearnedClause(out_learnt, out_btlevel);

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