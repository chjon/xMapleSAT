/***************************************************************************************[Solver.cc]
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

#include <math.h>

#include "core/Solver.h"
#include "mtl/Sort.h"

using namespace Minisat;

//=================================================================================================
// Options:


static const char* _cat = "CORE";

static BoolOption   opt_luby_restart (_cat, "luby",   "Use the Luby restart sequence", true);
static IntOption    opt_restart_first(_cat, "rfirst", "The base restart interval", 100, IntRange(1, INT32_MAX));
static DoubleOption opt_restart_inc  (_cat, "rinc",   "Restart interval increase factor", 2, DoubleRange(1, false, HUGE_VAL, false));

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

Solver::Solver()
    // Member variables
    : ok(true)
    , asynch_interrupt(false)

    // Parameters
    , luby_restart(opt_luby_restart)
    , restart_first(opt_restart_first)
    , restart_inc(opt_restart_inc)
    , conflict_budget(-1)
    , verbosity(0)

    // Statistics
    , solves(0)
    , starts(0)
    , conflicts(0)
    , simpDB_assigns(-1)
    , simpDB_props(0)

    // Solver components
    , assignmentTrail          (*this)
    , propagationQueue         (*this)
    , unitPropagator           (*this)
    , branchingHeuristicManager(*this)
    , clauseDatabase           (*this)
    , conflictAnalyzer         (*this)
{}

Solver::~Solver() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// PROBLEM SPECIFICATION

Var Solver::newVar(bool sign, bool dvar) {
    int v = assignmentTrail  .newVar();
    clauseDatabase           .newVar(v);
    propagationQueue         .newVar(v);
    unitPropagator           .newVar(v);
    branchingHeuristicManager.newVar(v, sign, dvar);
    conflictAnalyzer         .newVar(v);
    return v;
}

bool Solver::addClause(vec<Lit>& ps) {
    assert(assignmentTrail.decisionLevel() == 0);

    // Only add clause if the solver is in a consistent state
    if (!ok) return false;

    // Check if clause is satisfied and remove false/duplicate literals
    sort(ps);
    Lit p; int i, j;
    for (i = j = 0, p = lit_Undef; i < ps.size(); i++) {
        if (assignmentTrail.value(ps[i]) == l_True || ps[i] == ~p)
            return true;
        else if (assignmentTrail.value(ps[i]) != l_False && ps[i] != p)
            ps[j++] = p = ps[i];
    }
    ps.shrink(i - j);

    // Mark solver as inconsistent if the clause is empty
    if (ps.size() == 0) return ok = false;

    branchingHeuristicManager.handleEventInputClause(ps);

    // Add variable directly to trail if the clause is unit
    if (ps.size() == 1) {
        propagationQueue.enqueue(ps[0]);
        return ok = (unitPropagator.propagate() == CRef_Undef);
    }

    // Add the clause to the clause database
    clauseDatabase.addInputClause(ps);

    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

bool Solver::simplify(void) {
    assert(assignmentTrail.decisionLevel() == 0);

    // Don't need to simplify if the formula is already UNSAT
    if (!ok || unitPropagator.propagate() != CRef_Undef)
        return ok = false;

    if (assignmentTrail.nAssigns() == simpDB_assigns || (simpDB_props > 0))
        return true;

    // Remove satisfied clauses:
    clauseDatabase.removeSatisfied();
    
    // Add variables back to queue of decision variables
    branchingHeuristicManager.rebuildPriorityQueue();

    // Update stats (shouldn't depend on stats really, but it will do for now)
    simpDB_assigns = assignmentTrail.nAssigns();
    simpDB_props   = clauseDatabase.clauses_literals + clauseDatabase.learnts_literals;

    return true;
}

lbool Solver::search(int nof_conflicts) {
    assert(ok);
    int backtrack_level;
    int conflictC = 0;
    starts++;

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

            // Add the learnt clause to the clause database
            CRef cr = clauseDatabase.addLearntClause(learnt_clause);

            // First UIP learnt clauses are asserting after backjumping -- propagate!
            propagationQueue.enqueue(learnt_clause[0], cr);

        } else {
            // NO CONFLICT
            if (nof_conflicts >= 0 && conflictC >= nof_conflicts || !withinBudget()) {
                // Reached bound on number of conflicts:
                assignmentTrail.cancelUntil(0);
                return l_Undef;
            }

            // Simplify the set of problem clauses:
            if (assignmentTrail.decisionLevel() == 0 && !simplify())
                return l_False;

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
lbool Solver::solve_() {
    model.clear();
    conflict.clear();
    if (!ok) return l_False;

    solves++;

    clauseDatabase.init();

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
