/*************************************************************************************[ERSolver.cc]
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

#include "er/ERSolver.h"

using namespace Minisat;

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

ERSolver::ERSolver()
    : Solver()
    , erManager(*this)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

lbool ERSolver::search(int nof_conflicts) {
    assert(ok);
    int backtrack_level;
    int conflictC = 0;
    starts++;

#if ER_USER_GEN_LOCATION == ER_GEN_LOCATION_AFTER_RESTART
    // Generate extension variable definitions
    // Only try generating more extension variables if there aren't any buffered already
    erManager.checkGenerateDefinitions(conflicts);
#endif

#if ER_USER_ADD_LOCATION == ER_ADD_LOCATION_AFTER_RESTART
    // Add extension variables
    erManager.introduceExtVars();
#endif

    for (;;) {
        CRef confl = unitPropagator.propagate();
        branchingHeuristicManager.handleEventPropagated(conflicts, confl == CRef_Undef);

        if (confl != CRef_Undef) {
            // CONFLICT
            conflicts++; conflictC++;
            branchingHeuristicManager.handleEventConflicted(conflicts);

            // Check for root-level conflict
            if (assignmentTrail.decisionLevel() == 0) return l_False;

            // Generate a learnt clause from the conflict graph
            learnt_clause.clear();
            conflictAnalyzer.analyze(confl, learnt_clause, backtrack_level);

            // Backjump
            assignmentTrail.cancelUntil(backtrack_level);

            // EXTENDED RESOLUTION - substitute disjunctions with extension variables. This must be
            // called after backtracking because extension variables might need to be propagated.
            erManager.substitute(learnt_clause);

            // Add the learnt clause to the clause database
            CRef cr = clauseDatabase.addLearntClause(learnt_clause);

            // First UIP learnt clauses are asserting after backjumping -- propagate!
            propagationQueue.enqueue(learnt_clause[0], cr);

            // Filter clauses for clause selection
            if (cr != CRef_Undef)
                erManager.filterIncremental(cr);

        } else {
            // NO CONFLICT
            if (nof_conflicts >= 0 && conflictC >= nof_conflicts || !withinBudget()) {
                // Reached bound on number of conflicts:
                assignmentTrail.cancelUntil(0);
                return l_Undef;
            }

            // Simplify the set of problem clauses:
            if (assignmentTrail.decisionLevel() == 0) {
                if (!simplify()) return l_False;

            #if ER_USER_DELETE_HEURISTIC != ER_DELETE_HEURISTIC_NONE
                erManager.checkDeleteExtVars(conflicts);
            #endif
            }

            // Reduce the set of learnt clauses:
            clauseDatabase.checkReduceDB();

            Lit next = lit_Undef;
            while (assignmentTrail.decisionLevel() < assumptions.size()) {
                // Perform user provided assumption:
                Lit p = assumptions[assignmentTrail.decisionLevel()];
                if (assignmentTrail.value(p) == l_True) {
                    // Dummy decision level:
                    assignmentTrail.newDecisionLevel();
                } else if (assignmentTrail.value(p) == l_False) {
                    conflictAnalyzer.analyzeFinal(~p, conflict);
                    return l_False;
                } else {
                    next = p;
                    break;
                }
            }

            if (next == lit_Undef) {
                // New variable decision:
                next = branchingHeuristicManager.pickBranchLit();

                if (next == lit_Undef)
                    // Model found:
                    return l_True;

                // Update stats
                if (erManager.isExtVar(var(next)))
                    erManager.branchOnExt++;
            }

            // Increase decision level and enqueue 'next'
            assignmentTrail.newDecisionLevel();
            propagationQueue.enqueue(next);
        }
    }
}

/*
  Finite subsequences of the Luby-sequence:

  0: 1
  1: 1 1 2
  2: 1 1 2 1 1 2 4
  3: 1 1 2 1 1 2 4 1 1 2 1 1 2 4 8
  ...

 */

static double luby(double y, int x) {
    // Find the finite subsequence that contains index 'x', and the
    // size of that subsequence:
    int size, seq;
    for (size = 1, seq = 0; size < x+1; seq++, size = 2*size+1);

    while (size-1 != x){
        size = (size-1)>>1;
        seq--;
        x = x % size;
    }

    return pow(y, seq);
}

// NOTE: assumptions passed in member-variable 'assumptions'.
lbool ERSolver::solve_() {
    model.clear();
    conflict.clear();
    if (!ok) return l_False;

    solves++;

    // Initialize solver components
    clauseDatabase.init();
    erManager.init();

    lbool status = l_Undef;

    if (verbosity >= 1){
        printf("LBD Based Clause Deletion : %d\n", LBD_BASED_CLAUSE_DELETION);
        printf("Rapid Deletion : %d\n", RAPID_DELETION);
        printf("Almost Conflict : %d\n", ALMOST_CONFLICT);
        printf("Anti Exploration : %d\n", ANTI_EXPLORATION);
        printf("============================[ Search Statistics ]==============================\n");
        printf("| Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
        printf("|           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
        printf("===============================================================================\n");
    }

    // Search:
    int curr_restarts = 0;
    while (status == l_Undef){
        double rest_base = luby_restart ? luby(restart_inc, curr_restarts) : pow(restart_inc, curr_restarts);
        status = search(rest_base * restart_first);
        if (!withinBudget()) break;
        curr_restarts++;
    }

    if (verbosity >= 1)
        printf("===============================================================================\n");

    if (status == l_True) {
        // Extend & copy model:
        model.growTo(assignmentTrail.nVars());
        for (int i = 0; i < assignmentTrail.nVars(); i++) model[i] = assignmentTrail.value(i);
    } else if (status == l_False && conflict.size() == 0) {
        ok = false;
    }

    assignmentTrail.cancelUntil(0);
    return status;
}
