/*****************************************************************************[ConflictAnalyzer.cc]
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

#include "core/ConflictAnalyzer.h"
#include "core/Solver.h"

using namespace Minisat;

///////////////////////////////////////////////////////////////////////////////////////////////////
// OPTIONS

static const char* _cat = "CORE";

static IntOption opt_ccmin_mode(_cat, "ccmin-mode",  "Controls conflict clause minimization (0=none, 1=basic, 2=deep)", 2, IntRange(0, 2));

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

ConflictAnalyzer::ConflictAnalyzer(Solver& s)
    ////////////////////
    // Solver references
    : assignmentTrail(s.assignmentTrail)
    , branchingHeuristicManager(s.branchingHeuristicManager)
    , ca(s.ca)
    , solver(s)

    /////////////
    // Parameters
    , ccmin_mode(static_cast<ConflictClauseMinimizationMode>(static_cast<int>(opt_ccmin_mode)))

    //////////////////////
    // Temporary variables
    , counter(0)

    /////////////
    // Statistics
    , max_literals(0)
    , tot_literals(0)
{}


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
void ConflictAnalyzer::analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel, int& out_lbd) {
    // Generate conflict clause:
    firstUIPClause.clear();
    getFirstUIPClause(confl, firstUIPClause);
    max_literals += firstUIPClause.size();

    // Simplify conflict clause:
    simplifyClause(out_learnt, firstUIPClause);
    tot_literals += out_learnt.size();

    // Compute output LBD
    out_lbd = assignmentTrail.computeLBD(out_learnt);
    if (out_lbd <= 6 && out_learnt.size() <= 30) // Try further minimization?
        if (binResMinimize(out_learnt))
            out_lbd = assignmentTrail.computeLBD(out_learnt); // Recompute LBD if minimized.

    // Enforce watcher invariant
    enforceWatcherInvariant(out_learnt);

    // Find correct backtrack level
    out_btlevel = (out_learnt.size() == 1) ? 0 : assignmentTrail.level(var(out_learnt[1]));

    // Update data structures for branching heuristics
    branchingHeuristicManager.handleEventLearnedClause(out_learnt, seen, out_btlevel);

    // Clean up
    // TODO: can this be moved before updating the branching heuristic?
    // it is currently polluted by simplifyClause
    for (int j = 0; j < toClear.size(); j++)
        seen[toClear[j]] = false;
    toClear.clear();

    // Clear 'seen[]'
    for (int j = 0; j < firstUIPClause.size(); j++)
        seen[var(firstUIPClause[j])] = false;
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

    // Return an empty set of assumptions
    if (assignmentTrail.decisionLevel() == 0) return;

    // Mark the conflicting literal as seen
    seen[var(p)] = true;

    for (int i = assignmentTrail.nAssigns() - 1; i >= assignmentTrail.indexOfDecisionLevel(1); i--){
        Var x = var(assignmentTrail[i]);
        if (!seen[x]) continue;
        seen[x] = 0;

        if (assignmentTrail.reason(x) == CRef_Undef){
            assert(assignmentTrail.level(x) > 0);
            out_conflict.push(~assignmentTrail[i]);
        } else {
            Clause& c = ca[assignmentTrail.reason(x)];
            for (int j = c.size() == 2 ? 0 : 1; j < c.size(); j++)
                if (assignmentTrail.level(var(c[j])) > 0)
                    seen[var(c[j])] = 1;
        }
    }

    seen[var(p)] = 0;
}

template <class C>
static inline void enforceBinaryClauseInvariant(const AssignmentTrail& at, C& c) {
    if (c.size() == 2 && at.value(c[0]) == l_False) {
        // Special case for binary clauses - the first one has to be SAT
      assert(at.value(c[1]) != l_False);
      std::swap(c[0], c[1]);
    }
}

template <class C, class F>
static inline void forNonRootVariables(const AssignmentTrail& at, const C& c, int start, F f) {
    for (int i = start; i < c.size(); i++) {
        if (at.level(var(c[i])) != 0) {
            f(var(c[i]));
        }
    }
}

// pathCs[k] is the number of variables assigned at level k,
// it is initialized to 0 at the begining and reset to 0 after the function execution
bool ConflictAnalyzer::collectFirstUIP(CRef confl) {
    int totalPathCount = 0;

    // Mark distances of variables at the conflict
    forNonRootVariables(assignmentTrail, ca[confl], 0, [&](Var v) {
        // Mark distance of variable from the conflict
        var_iLevel_tmp[v] = 1;

        // Update data structures for exploring the conflict graph
        seen[v] = 1;
        pathCs[assignmentTrail.level(v)]++;
        totalPathCount++;
    });

    // Find every literal that participates in the conflict graph
    int max_distance = 1;
    involved_vars.clear();
    for (int i = assignmentTrail.nAssigns() - 1; totalPathCount > 0; i--) {
        Lit p = assignmentTrail[i];
        Var v = var(p);
        if (!seen[v]) continue;
        seen[v] = 0;

        // Keep track of variables that participate in the conflict graph
        involved_vars.push(var(p));

        // Don't explore backward past the first UIP for the level
        totalPathCount--;
        if (--pathCs[assignmentTrail.level(v)] == 0) continue;

        // Find the maximum path distance of a literal from the conflict
        const int reasonDist = var_iLevel_tmp[v] + 1;
        max_distance = std::max(max_distance, reasonDist);

        Clause& rc = ca[assignmentTrail.reason(v)];
        enforceBinaryClauseInvariant(assignmentTrail, rc);

        // Iterate over the non-root literals in the reason clause
        forNonRootVariables(assignmentTrail, rc, 1, [&](Var x) {
            // Mark distance of variable from the conflict
            var_iLevel_tmp[x] = seen[x] ? std::max(var_iLevel_tmp[x], reasonDist) : reasonDist;

            // Update data structures for exploring the conflict graph
            if (!seen[x]) {
                seen[x] = 1;
                totalPathCount++;
                pathCs[assignmentTrail.level(x)]++;
            }
        });
    }

    // Update activity_distance
    branchingHeuristicManager.updateActivityDistance(involved_vars, var_iLevel_tmp, max_distance);
    return true;
}

void ConflictAnalyzer::simpleAnalyze(
    CRef confl,
    vec<Lit>& out_learnt,
    vec<CRef>& reason_clause,
    bool True_confl,
    int trailRecord
) {
    int pathC = 0;
    Lit p = lit_Undef;
    int index = assignmentTrail.nAssigns() - 1;

    do {
        if (confl != CRef_Undef){
            reason_clause.push(confl);
            Clause& c = ca[confl];
	    if (p != lit_Undef or True_confl == true) 
	      enforceBinaryClauseInvariant(assignmentTrail, c);

            // if True_confl==true, then choose p begin with the 1th index of c;
            for (int j = (p == lit_Undef && True_confl == false) ? 0 : 1; j < c.size(); j++){
                Lit q = c[j];
                if (!seen[var(q)]){
                    seen[var(q)] = 1;
                    pathC++;
                }
            }
        }
        else {
            out_learnt.push(~p);
        }
        // if not break, while() will come to the index of trail blow 0, and fatal error occur;
        if (pathC == 0) break;
        // Select next clause to look at:
        while (!seen[var(assignmentTrail[index--])]);
        // if the reason cr from the 0-level assigned var, we must break avoid move forth further;
        // but attention that maybe seen[x]=1 and never be clear. However makes no matter;
        if (trailRecord > index + 1) break;
        p = assignmentTrail[index + 1];
        confl = assignmentTrail.reason(var(p));
        seen[var(p)] = 0;
        pathC--;

    } while (pathC >= 0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

// Try further learnt clause minimization by means of binary clause resolution.
bool ConflictAnalyzer::binResMinimize(vec<Lit>& out_learnt) {
    // Preparation: remember which false variables we have in 'out_learnt'.
    counter++;
    for (int i = 1; i < out_learnt.size(); i++)
        seen2[var(out_learnt[i])] = counter;

    // Get the list of binary clauses containing 'out_learnt[0]'.
    const vec<Watcher>& ws = solver.unitPropagator.getBinWatchers(~out_learnt[0]);

    int to_remove = 0;
    for (int i = 0; i < ws.size(); i++) {
        Lit the_other = ws[i].blocker;
        // Does 'the_other' appear negatively in 'out_learnt'?
        if (seen2[var(the_other)] == counter && assignmentTrail.value(the_other) == l_True){
            to_remove++;
            seen2[var(the_other)] = counter - 1; // Remember to remove this variable.
        }
    }

    // Shrink.
    if (to_remove > 0) {
        int last = out_learnt.size() - 1;
        for (int i = 1; i < out_learnt.size() - to_remove; i++)
            if (seen2[var(out_learnt[i])] != counter)
                out_learnt[i--] = out_learnt[last--];
        out_learnt.shrink(to_remove);
    }
    return to_remove != 0;
}

// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool ConflictAnalyzer::litRedundant(Lit p, uint32_t abstract_levels) {
    // Initialize local data structures
    const int top = toClear.size();
    workStack.push(assignmentTrail.reason(var(p)));

    // Iterate through reason clauses
    while (workStack.size() > 0){
        assert(workStack.last() != CRef_Undef);
        Clause& c = ca[workStack.last()]; workStack.pop();

        // Special handling for binary clauses like in 'analyze()'.
        if (c.size() == 2 && assignmentTrail.value(c[0]) == l_False){
            assert(assignmentTrail.value(c[1]) == l_True);
            std::swap(c[0], c[1]);
        }

        // Iterate through unique reason variables that are not assigned at the root level
        for (int i = 1; i < c.size(); i++) {
            const Var v = var(c[i]);
            if (seen[v] || assignmentTrail.level(v) == 0) continue;

            // Clean up and abort (don't remove literal) if:
            //     1. a decision variable is the reason variable OR
            //     2. the literal is at a level that cannot be removed later
            const CRef reason = assignmentTrail.reason(v);
            if (reason == CRef_Undef ||
                (assignmentTrail.abstractLevel(v) & abstract_levels) == 0
            ) {
                // Clean up
                for (int j = top; j < toClear.size(); j++)
                    seen[toClear[j]] = false;
                toClear.shrink(toClear.size() - top);
                workStack.clear();

                return false;
            }

            // Mark variable as seen and add its reason clause to the work queue
            seen[v] = true;
            toClear.push(v);
            workStack.push(reason);
        }
    }

    return true;
}

inline void ConflictAnalyzer::getFirstUIPClause(CRef confl, vec<Lit>& out_learnt) {
    Lit p = lit_Undef;
    int pathC = 0;
    int index = assignmentTrail.nAssigns() - 1;
    out_learnt.push(); // (leave room for the asserting literal)

    do {
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

        // For binary clauses, we don't rearrange literals in propagate(), so check and make sure the first is an implied lit.
        if (p != lit_Undef && c.size() == 2 && assignmentTrail.value(c[0]) == l_False) {
            assert(assignmentTrail.value(c[1]) == l_True);
            std::swap(c[0], c[1]);
        }

        solver.clauseDatabase.handleEventClauseInConflictGraph(confl, solver.conflicts);

        // Iterate through every literal that participates in the conflict graph
        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++) {
            Lit q = c[j];
            if (seen[var(q)] || assignmentTrail.level(var(q)) == 0) continue;

            // Mark variable as seen
            branchingHeuristicManager.handleEventLitInConflictGraph(q, solver.conflicts);
            seen[var(q)] = 1;

            // Increment the number of paths if the variable is assigned at the decision level that resulted in conflict
            if (assignmentTrail.level(var(q)) >= assignmentTrail.decisionLevel())
                pathC++;
            // Add literals from earlier decision levels to the conflict clause
            else
                out_learnt.push(q);
        }
        
        // Select next clause to look at:
        while (!seen[var(assignmentTrail[index--])]);
        p     = assignmentTrail[index+1];
        confl = assignmentTrail.reason(var(p));

        // Mark variable as unseen: it is either at or after the first UIP
        seen[var(p)] = 0;
        pathC--;

    } while (pathC > 0);

    // Add first UIP literal at index 0
    out_learnt[0] = ~p;
    seen[var(p)] = true;

    // Note: at this point, seen[v] is true iff v is in the learnt clause
}

inline void ConflictAnalyzer::simplifyClauseDeep(vec<Lit>& simplified, const vec<Lit>& toSimplify) {
    // Initialize abstraction of levels involved in conflict
    uint32_t abstract_level = 0;
    for (int i = 1; i < toSimplify.size(); i++)
        abstract_level |= assignmentTrail.abstractLevel(var(toSimplify[i]));

    // Copy non-redundant literals
    simplified.push(toSimplify[0]);
    for (int i = 1; i < toSimplify.size(); i++) {
        if (// Keep decision literals
            assignmentTrail.reason(var(toSimplify[i])) == CRef_Undef ||
            
            // Keep literals that are not redundant
            !litRedundant(toSimplify[i], abstract_level)
        ) simplified.push(toSimplify[i]);
    }
}

inline bool ConflictAnalyzer::reasonSubsumed(const Clause& c) {
    // Iterate through every variable in the reason clause, ignoring the propagated variable
    for (int k = c.size() == 2 ? 0 : 1; k < c.size(); k++) {
        // If a non-root variable is not in the learnt clause, the reason clause is not subsumed!
        if (!seen[var(c[k])] && assignmentTrail.level(var(c[k])) > 0)
            return false;
    }

    return true;
}

inline void ConflictAnalyzer::simplifyClauseBasic(vec<Lit>& simplified, const vec<Lit>& toSimplify) {
    // Iterate through every variable in the learnt clause (excluding the asserting literal)
    simplified.push(toSimplify[0]);
    for (int i = 1; i < toSimplify.size(); i++){
        const CRef reason = assignmentTrail.reason(var(toSimplify[i]));
        if (// Keep decision variables
            reason == CRef_Undef ||

            // Keep variables whose reason clauses are not subsumed by the learnt clause
            !reasonSubsumed(ca[reason])
        ) simplified.push(toSimplify[i]);
    }
}

inline void ConflictAnalyzer::simplifyClause(vec<Lit>& simplified, const vec<Lit>& toSimplify) {
    switch (ccmin_mode) {
        case ConflictClauseMinimizationMode::DEEP:  return simplifyClauseDeep(simplified, toSimplify);
        case ConflictClauseMinimizationMode::BASIC: return simplifyClauseBasic(simplified, toSimplify);
        default: return toSimplify.copyTo(simplified);
    }
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
