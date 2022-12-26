/***************************************************************************************[Solver.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

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

#if ! LBD_BASED_CLAUSE_DELETION
static DoubleOption opt_clause_decay (_cat, "cla-decay", "The clause activity decay factor", 0.999, DoubleRange(0, false, 1, false));
#endif
static BoolOption   opt_luby_restart (_cat, "luby",   "Use the Luby restart sequence", true);
static IntOption    opt_restart_first(_cat, "rfirst", "The base restart interval", 100, IntRange(1, INT32_MAX));
static DoubleOption opt_restart_inc  (_cat, "rinc",   "Restart interval increase factor", 2, DoubleRange(1, false, HUGE_VAL, false));


//=================================================================================================
// Constructor/Destructor:


Solver::Solver() :

    // Parameters (user settable):
    //
    verbosity        (0)
#if ! LBD_BASED_CLAUSE_DELETION
  , clause_decay     (opt_clause_decay)
#endif
  , luby_restart     (opt_luby_restart)
  , restart_first    (opt_restart_first)
  , restart_inc      (opt_restart_inc)

    // Parameters (the rest):
    //
  , learntsize_factor((double)1/(double)3), learntsize_inc(1.1)

    // Parameters (experimental):
    //
  , learntsize_adjust_start_confl (100)
  , learntsize_adjust_inc         (1.5)

    // Statistics: (formerly in 'SolverStats')
    //
  , solves(0), starts(0), conflicts(0)
  , lbd_calls(0)

  , ok                 (true)
#if ! LBD_BASED_CLAUSE_DELETION
  , cla_inc            (1)
#endif
  , simpDB_assigns     (-1)
  , simpDB_props       (0)
  , progress_estimate  (0)

    // Resource constraints:
    //
  , conflict_budget    (-1)
  , propagation_budget (-1)
  , asynch_interrupt   (false)

  // Solver components
  , assignmentTrail(this)
  , propagationQueue(this)
  , unitPropagator(this)
  , branchingHeuristicManager(this)
  , clauseDatabase(this)
  , conflictAnalyzer(this)
{}


Solver::~Solver()
{
}


//=================================================================================================
// Minor methods:


// Creates a new SAT variable in the solver. If 'decision' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool sign, bool dvar) {
    int v = variableDatabase .newVar();
    assignmentTrail          .newVar(v);
    propagationQueue         .newVar(v);
    unitPropagator           .newVar(v);
    branchingHeuristicManager.newVar(v, sign, dvar);
    conflictAnalyzer         .newVar(v);
    seen.push(0);
    lbd_seen.push(0);
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
        if (variableDatabase.value(ps[i]) == l_True || ps[i] == ~p)
            return true;
        else if (variableDatabase.value(ps[i]) != l_False && ps[i] != p)
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

int min(int a, int b) {
    return a < b ? a : b;
}

/*_________________________________________________________________________________________________
|
|  simplify : [void]  ->  [bool]
|  
|  Description:
|    Simplify the clause database according to the current top-level assigment. Currently, the only
|    thing done here is the removal of satisfied clauses, but more things can be put here.
|________________________________________________________________________________________________@*/
bool Solver::simplify() {
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

    simpDB_assigns = assignmentTrail.nAssigns();
    simpDB_props   = clauseDatabase.clauses_literals + clauseDatabase.learnts_literals;   // (shouldn't depend on stats really, but it will do for now)

    return true;
}

/*_________________________________________________________________________________________________
|
|  search : (nof_conflicts : int) (params : const SearchParams&)  ->  [lbool]
|  
|  Description:
|    Search for a model the specified number of conflicts. 
|    NOTE! Use negative value for 'nof_conflicts' indicate infinity.
|  
|  Output:
|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
|________________________________________________________________________________________________@*/
lbool Solver::search(int nof_conflicts) {
    assert(ok);
    int         backtrack_level;
    int         conflictC = 0;
    vec<Lit>    learnt_clause;
    starts++;

    for (;;) {
        CRef confl = unitPropagator.propagate();

        branchingHeuristicManager.handleEventPropagated(conflicts, confl == CRef_Undef);

        if (confl != CRef_Undef){
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

            // Update clause activity
            // TODO: the solver appears to perform better if we move this block after enqueuing
            if (cr != CRef_Undef) {
            #if LBD_BASED_CLAUSE_DELETION
                Clause& clause = ca[cr];
                clause.activity() = lbd(clause);
            #else
                claBumpActivity(ca[cr]);
            #endif
            }

            // First UIP learnt clauses are asserting after backjumping -- propagate!
            propagationQueue.enqueue(learnt_clause[0], cr);

#if BRANCHING_HEURISTIC == VSIDS
            branchingHeuristicManager.decayActivityVSIDS();
#endif
#if ! LBD_BASED_CLAUSE_DELETION
            claDecayActivity();
#endif

            if (--learntsize_adjust_cnt == 0){
                learntsize_adjust_confl *= learntsize_adjust_inc;
                learntsize_adjust_cnt    = (int)learntsize_adjust_confl;
#if ! RAPID_DELETION
                max_learnts             *= learntsize_inc;
#endif

                if (verbosity >= 1)
                    printf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %6.3f %% |\n", 
                           (int)conflicts, 
                           nFreeVars(), clauseDatabase.nClauses(), (int)clauseDatabase.clauses_literals, 
                           (int)max_learnts, clauseDatabase.nLearnts(), (double)clauseDatabase.learnts_literals/clauseDatabase.nLearnts(), assignmentTrail.progressEstimate()*100);
            }

        } else {
            // NO CONFLICT
            if (nof_conflicts >= 0 && conflictC >= nof_conflicts || !withinBudget()) {
                // Reached bound on number of conflicts:
                progress_estimate = assignmentTrail.progressEstimate();
                assignmentTrail.cancelUntil(0);
                return l_Undef;
            }

            // Simplify the set of problem clauses:
            if (assignmentTrail.decisionLevel() == 0 && !simplify())
                return l_False;

            // Reduce the set of learnt clauses:
            if (clauseDatabase.nLearnts() - assignmentTrail.nAssigns() >= max_learnts) {
                clauseDatabase.reduceDB();
#if RAPID_DELETION
                max_learnts += 500;
#endif
            }

            Lit next = lit_Undef;
            while (assignmentTrail.decisionLevel() < assumptions.size()) {
                // Perform user provided assumption:
                Lit p = assumptions[assignmentTrail.decisionLevel()];
                if (variableDatabase.value(p) == l_True) {
                    // Dummy decision level:
                    assignmentTrail.newDecisionLevel();
                } else if (variableDatabase.value(p) == l_False) {
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

static double luby(double y, int x){

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

#if RAPID_DELETION
    max_learnts               = 2000;
#else
    max_learnts               = nClauses() * learntsize_factor;
#endif
    learntsize_adjust_confl   = learntsize_adjust_start_confl;
    learntsize_adjust_cnt     = (int)learntsize_adjust_confl;
    lbool   status            = l_Undef;

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
        model.growTo(variableDatabase.nVars());
        for (int i = 0; i < variableDatabase.nVars(); i++) model[i] = variableDatabase.value(i);
    } else if (status == l_False && conflict.size() == 0) {
        ok = false;
    }

    assignmentTrail.cancelUntil(0);
    return status;
}
