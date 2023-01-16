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

lbool ERSolver::search(int& nof_conflicts) {
    assert(ok);
    int         backtrack_level;
    int         lbd;
    vec<Lit>    learnt_clause;
    starts++;

    restartHeuristicManager.handleEventRestarted(nof_conflicts);

    #if ER_USER_GEN_LOCATION == ER_GEN_LOCATION_AFTER_RESTART
        // Generate extension variable definitions
        // Only try generating more extension variables if there aren't any buffered already
        erManager.checkGenerateDefinitions(conflicts);
    #endif

    #if ER_USER_ADD_LOCATION == ER_ADD_LOCATION_AFTER_RESTART
        // Add extension variables
        erManager.introduceExtVars();
    #endif

    // Simplify
    if (conflicts >= static_cast<unsigned int>(curSimplify * nbconfbeforesimplify)) {
        if (!simplifyAll()) return l_False;
        curSimplify = (conflicts / nbconfbeforesimplify) + 1;
        nbconfbeforesimplify += incSimplify;
    }

    for (;;){
        CRef confl = unitPropagator.propagate();

        if (confl != CRef_Undef) {
            // CONFLICT
            conflicts++;
            if (assignmentTrail.decisionLevel() == 0) return l_False;

            clauseDatabase.handleEventConflicted(conflicts);
            branchingHeuristicManager.handleEventConflicted(confl, conflicts);

            learnt_clause.clear();
            conflictAnalyzer.analyze(confl, learnt_clause, backtrack_level, lbd);
            assignmentTrail.cancelUntilLevel(backtrack_level);

            // EXTENDED RESOLUTION - substitute disjunctions with extension variables. This must be
            // called after backtracking because extension variables might need to be propagated.
            erManager.substitute(learnt_clause);

            lbd--;
            restartHeuristicManager.handleEventLearntClause(lbd);

            CRef cr = clauseDatabase.addLearntClause(learnt_clause, lbd, conflicts);
            propagationQueue.enqueue(learnt_clause[0], cr);

            // Filter clauses for clause selection
            if (cr != CRef_Undef)
                erManager.filterIncremental(cr);

        #if ER_USER_GEN_LOCATION == ER_GEN_LOCATION_AFTER_CONFLICT
            // Generate extension variable definitions
            // Only try generating more extension variables if there aren't any buffered already
            erManager.checkGenerateDefinitions(conflicts);
        #endif

        #if ER_USER_ADD_LOCATION == ER_ADD_LOCATION_AFTER_CONFLICT
            // Add extension variables
            erManager.introduceExtVars();
        #endif

        #if ER_ENABLE_GLUCOSER
            erManager.generateLER();
            erManager.introduceExtVars(ERManager::HeuristicType::LER);
        #endif

            // Output to proof file
            proofLogger.addClause(learnt_clause);
        } else {
            // NO CONFLICT
            if (restartHeuristicManager.shouldRestart() || !withinBudget()) {
                // Reached bound on number of conflicts:
                assignmentTrail.cancelUntilLevel(0);
                nof_conflicts = restartHeuristicManager.getConflictBudget();
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
            clauseDatabase.checkReduceDB(conflicts);

            // New variable decision:
            Lit next = branchingHeuristicManager.pickBranchLit();
            if (next == lit_Undef)
                // Model found:
                return l_True;

            // Update stats
            if (erManager.isExtVar(var(next)))
                erManager.branchOnExt++;

            // Increase decision level and enqueue 'next'
            assignmentTrail.newDecisionLevel();
            propagationQueue.enqueue(next);
        }
    }
}

// NOTE: assumptions passed in member-variable 'assumptions'.
lbool ERSolver::solve_() {
    model.clear();
    conflict.clear();
    if (!ok) return l_False;

    solves++;

    // Initialize solver components
    erManager.init();

    lbool status = l_Undef;

    if (verbosity >= 1){
        printf("c ============================[ Search Statistics ]==============================\n");
        printf("c | Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
        printf("c |           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
        printf("c ===============================================================================\n");
    }

    branchingHeuristicManager.VSIDS = true;
    int init = 10000;
    while (status == l_Undef && init > 0 && withinBudget())
        status = search(init);
    branchingHeuristicManager.VSIDS = false;

    // Search:
    while (status == l_Undef && withinBudget()) {
        // Periodically switch branching heuristic
        branchingHeuristicManager.checkSwitchHeuristic(unitPropagator.propagations);

        // Compute the next number of conflicts before restart
        int numConflictsBeforeRestart = restartHeuristicManager.getRestartConflicts();

        // Search
        status = search(numConflictsBeforeRestart);
    }

    if (verbosity >= 1)
        printf("c ===============================================================================\n");

    if (status == l_False) proofLogger.flush();

    if (status == l_True) {
        // Extend & copy model:
        model.growTo(assignmentTrail.nVars());
        for (int i = 0; i < assignmentTrail.nVars(); i++) model[i] = assignmentTrail.value(i);
    } else if (status == l_False && conflict.size() == 0) {
        ok = false;
    }

    assignmentTrail.cancelUntilLevel(0);
    return status;
}