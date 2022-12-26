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
static DoubleOption  opt_clause_decay      (_cat, "cla-decay",   "The clause activity decay factor",              0.999,    DoubleRange(0, false, 1, false));
#endif
static DoubleOption  opt_random_var_freq   (_cat, "rnd-freq",    "The frequency with which the decision heuristic tries to choose a random variable", 0, DoubleRange(0, true, 1, true));
static IntOption     opt_ccmin_mode        (_cat, "ccmin-mode",  "Controls conflict clause minimization (0=none, 1=basic, 2=deep)", 2, IntRange(0, 2));
static BoolOption    opt_rnd_init_act      (_cat, "rnd-init",    "Randomize the initial activity", false);
static BoolOption    opt_luby_restart      (_cat, "luby",        "Use the Luby restart sequence", true);
static IntOption     opt_restart_first     (_cat, "rfirst",      "The base restart interval", 100, IntRange(1, INT32_MAX));
static DoubleOption  opt_restart_inc       (_cat, "rinc",        "Restart interval increase factor", 2, DoubleRange(1, false, HUGE_VAL, false));


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
  , ccmin_mode       (opt_ccmin_mode)
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
  , max_literals(0), tot_literals(0)

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
    seen     .push(0);
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
        assignmentTrail.uncheckedEnqueue(ps[0]);
        return ok = (unitPropagator.propagate() == CRef_Undef);
    }

    // Add the clause to the clause database
    clauseDatabase.addInputClause(ps);

    return true;
}

//=================================================================================================
// Major methods:

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
void Solver::analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel) {
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
            c.activity() = lbd(c);
#else
        if (c.learnt())
            claBumpActivity(c);
#endif

        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];

            if (!seen[var(q)] && assignmentTrail.level(var(q)) > 0){
                branchingHeuristicManager.handleEventLitInConflictGraph(q, conflicts);
                seen[var(q)] = 1;
                if (assignmentTrail.level(var(q)) >= assignmentTrail.decisionLevel())
                    pathC++;
                else
                    out_learnt.push(q);
            }
        }
        
        // Select next clause to look at:
        while (!seen[var(assignmentTrail[index--])]);
        p     = assignmentTrail[index+1];
        confl = assignmentTrail.reason(var(p));
        seen[var(p)] = 0;
        pathC--;

    } while (pathC > 0);
    out_learnt[0] = ~p;

    // Simplify conflict clause:
    //
    int i, j;
    out_learnt.copyTo(analyze_toclear);
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
                    if (!seen[var(c[k])] && assignmentTrail.level(var(c[k])) > 0){
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
    for (int j = 0; j < analyze_toclear.size(); j++)
        seen[var(analyze_toclear[j])] = 0;
}


// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, uint32_t abstract_levels) {
    analyze_stack.clear(); analyze_stack.push(p);
    int top = analyze_toclear.size();
    while (analyze_stack.size() > 0){
        assert(assignmentTrail.reason(var(analyze_stack.last())) != CRef_Undef);
        Clause& c = ca[assignmentTrail.reason(var(analyze_stack.last()))]; analyze_stack.pop();

        for (int i = 1; i < c.size(); i++){
            Lit p = c[i];
            if (seen[var(p)] || assignmentTrail.level(var(p)) == 0) continue;

            if (assignmentTrail.reason(var(p)) != CRef_Undef && (assignmentTrail.abstractLevel(var(p)) & abstract_levels) != 0){
                seen[var(p)] = 1;
                analyze_stack.push(p);
                analyze_toclear.push(p);
            } else {
                for (int j = top; j < analyze_toclear.size(); j++)
                    seen[var(analyze_toclear[j])] = 0;
                analyze_toclear.shrink(analyze_toclear.size() - top);
                return false;
            }
        }
    }

    return true;
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
void Solver::analyzeFinal(Lit p, vec<Lit>& out_conflict) {
    out_conflict.clear();
    out_conflict.push(p);

    if (assignmentTrail.decisionLevel() == 0)
        return;

    seen[var(p)] = 1;

    for (int i = assignmentTrail.nAssigns() - 1; i >= assignmentTrail.indexOfDecisionLevel(1); i--) {
        Var x = var(assignmentTrail[i]);
        if (seen[x]) {
            if (assignmentTrail.reason(x) == CRef_Undef){
                assert(assignmentTrail.level(x) > 0);
                out_conflict.push(~assignmentTrail[i]);
            } else {
                Clause& c = ca[assignmentTrail.reason(x)];
                for (int j = 1; j < c.size(); j++)
                    if (assignmentTrail.level(var(c[j])) > 0)
                        seen[var(c[j])] = 1;
            }
            seen[x] = 0;
        }
    }

    seen[var(p)] = 0;
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
            analyze(confl, learnt_clause, backtrack_level);

            // Backjump
            assignmentTrail.cancelUntil(backtrack_level);

            // Add the learnt clause to the clause database
            CRef cr = clauseDatabase.addLearntClause(learnt_clause);

            // Update clause activity
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
                    analyzeFinal(~p, conflict);
                    return l_False;
                } else {
                    next = p;
                    break;
                }
            }

            if (next == lit_Undef){
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


    if (status == l_True){
        // Extend & copy model:
        model.growTo(variableDatabase.nVars());
        for (int i = 0; i < variableDatabase.nVars(); i++) model[i] = variableDatabase.value(i);
    }else if (status == l_False && conflict.size() == 0)
        ok = false;

    assignmentTrail.cancelUntil(0);
    return status;
}

//=================================================================================================
// Garbage Collection methods:

void Solver::relocAll(ClauseAllocator& to) {
    // All watchers:
    unitPropagator.relocAll(to);

    // All reasons:
    assignmentTrail.relocAll(to);

    // All clauses:
    clauseDatabase.relocAll(to);
}
