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

static const char* _cat = "CORE";

static IntOption opt_ccmin_mode (_cat, "ccmin-mode", "Controls conflict clause minimization (0=none, 1=basic, 2=deep)", 2, IntRange(0, 2));

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

    /////////////
    // Statistics

    , max_literals(0)
    , tot_literals(0)

{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// PUBLIC API

void ConflictAnalyzer::analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel) {
    // Generate conflict clause:
    firstUIPClause.clear();
    getFirstUIPClause(confl, firstUIPClause);
    max_literals += firstUIPClause.size();

    // Simplify conflict clause:
    simplifyClause(out_learnt, firstUIPClause);
    tot_literals += out_learnt.size();

    // Enforce watcher invariant
    enforceWatcherInvariant(out_learnt);

    // Find backtrack level:
    out_btlevel = (out_learnt.size() == 1) ? 0 : assignmentTrail.level(var(out_learnt[1]));

    // Update data structures for branching heuristics
    branchingHeuristicManager.handleEventLearnedClause(out_learnt, seen);

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

void ConflictAnalyzer::analyzeFinal(Lit p, vec<Lit>& out_conflict) {
    out_conflict.clear();
    out_conflict.push(p);

    // Return an empty set of assumptions
    if (assignmentTrail.decisionLevel() == 0) return;

    // Mark the conflicting literal as seen
    seen[var(p)] = true;

    // Iterate through every variable in the conflict graph
    for (int i = assignmentTrail.nAssigns() - 1; i >= assignmentTrail.indexOfDecisionLevel(1); i--) {
        Var x = var(assignmentTrail[i]);
        if (!seen[x]) continue;

        if (assignmentTrail.reason(x) == CRef_Undef){
            // Add variable to conflict clause
            assert(assignmentTrail.level(x) > 0);
            out_conflict.push(~assignmentTrail[i]);
        } else {
            // Queue variables to be analyzed
            Clause& c = ca[assignmentTrail.reason(x)];
            for (int j = 1; j < c.size(); j++)
                if (assignmentTrail.level(var(c[j])) > 0)
                    seen[var(c[j])] = true;
        }
        
        seen[x] = false;
    }

    // Clear seen variables
    seen[var(p)] = false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

bool ConflictAnalyzer::litRedundant(Lit p, uint32_t abstract_levels) {
    // Initialize local data structures
    const int top = toClear.size();
    workStack.push(assignmentTrail.reason(var(p)));

    // Iterate through reason clauses
    while (workStack.size() > 0) {
        assert(workStack.last() != CRef_Undef);
        Clause& c = ca[workStack.last()]; workStack.pop();

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

    // TODO: bug? why not clean up the changes to seen?

    // Clean up
    workStack.clear();

    return true;
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
    for (int k = 1; k < c.size(); k++) {
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

inline void ConflictAnalyzer::getFirstUIPClause(CRef confl, vec<Lit>& learntClause) {
    // Initialize local data structures
    Lit p = lit_Undef;
    int pathC = 0;
    int index = assignmentTrail.nAssigns() - 1;
    learntClause.push(); // (leave room for the asserting literal)
    
    // Iterate through every clause that participates in the conflict graph, backtracking until the first UIP
    do {
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

        solver.clauseDatabase.handleEventClauseInConflictGraph(c);

        // Iterate through every literal that participates in the conflict graph
        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];
            if (seen[var(q)] || assignmentTrail.level(var(q)) == 0) continue;

            // Mark variable as seen
            branchingHeuristicManager.handleEventLitInConflictGraph(q, solver.conflicts);
            seen[var(q)] = true;

            // Increment the number of paths if the variable is assigned at the decision level that resulted in conflict
            if (assignmentTrail.level(var(q)) == assignmentTrail.decisionLevel())
                pathC++;
            // Add literals from earlier decision levels to the conflict clause
            else
                learntClause.push(q);
        }
        
        // Select next clause to look at:
        while (!seen[var(assignmentTrail[index--])]);
        p     = assignmentTrail[index+1];
        confl = assignmentTrail.reason(var(p));

        // Mark variable as unseen: it is either at or after the first UIP
        seen[var(p)] = false;
        pathC--;

    } while (pathC > 0);

    // Add first UIP literal at index 0
    learntClause[0] = ~p;
    seen[var(p)] = true;

    // Note: at this point, seen[v] is true iff v is in the learnt clause
}

inline void ConflictAnalyzer::simplifyClause(vec<Lit>& simplified, const vec<Lit>& toSimplify) {
    switch (ccmin_mode) {
        case ConflictClauseMinimizationMode::DEEP:  return simplifyClauseDeep(simplified, toSimplify);
        case ConflictClauseMinimizationMode::BASIC: return simplifyClauseBasic(simplified, toSimplify);
        default: return;
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