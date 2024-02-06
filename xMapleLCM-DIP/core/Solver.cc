/***************************************************************************************[Solver.cc]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010, Niklas Sorensson
 
Chanseok Oh's MiniSat Patch Series -- Copyright (c) 2015, Chanseok Oh
 
Maple_LCM, Based on MapleCOMSPS_DRUP -- Copyright (c) 2017, Mao Luo, Chu-Min LI, Fan Xiao: implementing a learnt clause minimisation approach
Reference: M. Luo, C.-M. Li, F. Xiao, F. Manya, and Z. L. , “An effective learnt clause minimization approach for cdcl sat solvers,” in IJCAI-2017, 2017, pp. to–appear.

Maple_LCM_Dist, Based on Maple_LCM -- Copyright (c) 2017, Fan Xiao, Chu-Min LI, Mao Luo: using a new branching heuristic called Distance at the beginning of search
 
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

#include <signal.h>
#include <unistd.h>

#include "mtl/Sort.h"
#include "core/Solver.h"

using namespace Minisat;
using namespace std;

static const char* _cat = "CORE";
static BoolOption   opt_produce_proof      (_cat, "produce-proof", "Produce DRAT proof in file proof.txt.", false);

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

Solver::Solver()
    // Member variables
    : ok(true)
    , asynch_interrupt(false)

    // Parameters
    , conflict_budget(-1)
    , curSimplify(1)
    , nbconfbeforesimplify(1000)
    , incSimplify(1000)
    , verbosity(0)
    , produce_proof(opt_produce_proof)
      
    // Statistics: (formerly in 'SolverStats')
    //
    , solves(0)
    , starts(0)
    , conflicts(0)
    , dip_conflicts(0)  
    , simpDB_assigns(-1)
    , simpDB_props(0)

    // Solver components
    , assignmentTrail          (*this)
    , branchingHeuristicManager(*this)
    , clauseDatabase           (*this)
    , conflictAnalyzer         (*this)
    , propagationQueue         (*this)
    , restartHeuristicManager  (*this)
    , unitPropagator           (*this)
{
  if (produce_proof)
    proofLogger.drup_file = fopen("proof.txt","w");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// SOLVING FUNCTIONS

lbool Solver::search(int& nof_conflicts) {
    assert(ok);
    int         backtrack_level;
    int         lbd;
    vec<Lit>    learnt_clause;
    starts++;

    restartHeuristicManager.handleEventRestarted(nof_conflicts);

    // Simplify
    if (conflicts >= static_cast<unsigned int>(curSimplify * nbconfbeforesimplify)) {
        if (!simplifyAll()) return l_False;
        curSimplify = (conflicts / nbconfbeforesimplify) + 1;
        nbconfbeforesimplify += incSimplify;
    }

    propagationQueue.updatePrioritizationMode(
        branchingHeuristicManager.DISTANCE,
        branchingHeuristicManager.getActivity()
    );

    for (;;){
        CRef confl = unitPropagator.propagate();

        if (confl != CRef_Undef) {
            // CONFLICT
            conflicts++;
            if (assignmentTrail.decisionLevel() == 0) return l_False;

            clauseDatabase.handleEventConflicted(conflicts);
            branchingHeuristicManager.handleEventConflicted(confl, conflicts);

            learnt_clause.clear();
	    vec<Lit> learnt_clause_2;
            conflictAnalyzer.analyze(confl, learnt_clause, backtrack_level, lbd, learnt_clause_2);
            assignmentTrail.cancelUntilLevel(backtrack_level);

            lbd--;
            restartHeuristicManager.handleEventLearntClause(lbd);

            CRef cr = clauseDatabase.addLearntClause(learnt_clause, lbd, conflicts);
            propagationQueue.enqueue(learnt_clause[0], cr);

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
            if (assignmentTrail.decisionLevel() == 0 && !simplify())
                return l_False;

            // Reduce the set of learnt clauses:
            clauseDatabase.checkReduceDB(conflicts);

            // New variable decision:
            Lit next = branchingHeuristicManager.pickBranchLit();
            if (next == lit_Undef)
                // Model found:
                return l_True;

            // Increase decision level and enqueue 'next'
            assignmentTrail.newDecisionLevel();
            propagationQueue.enqueue(next);
        }
    }
}

// NOTE: assumptions passed in member-variable 'assumptions'.
lbool Solver::solve_() {
    model.clear();
    conflict.clear();
    if (!ok) return l_False;

    solves++;

    lbool status = l_Undef;

    if (verbosity >= 1){
        printf("c ============================[ Search Statistics ]==============================\n");
        printf("c | Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
        printf("c |           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
        printf("c ===============================================================================\n");
    }

#if BRANCHING_HEURISTIC == BRANCHING_HEURISTIC_DYNAMIC
    branchingHeuristicManager.VSIDS = true;
#endif
    int init = 10000;
    while (status == l_Undef && init > 0 && withinBudget())
        status = search(init);
#if BRANCHING_HEURISTIC == BRANCHING_HEURISTIC_DYNAMIC
    branchingHeuristicManager.VSIDS = false;
#endif

    // Search:
    while (status == l_Undef && withinBudget()) {
#if BRANCHING_HEURISTIC == BRANCHING_HEURISTIC_DYNAMIC
        // Periodically switch branching heuristic
        branchingHeuristicManager.checkSwitchHeuristic(unitPropagator.propagations);
#endif

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

// simplify All
//

void Solver::simplifyLearnt(Clause& c) {
    const int trailRecord = assignmentTrail.nAssigns(); // record the start pointer

    vec<Lit> falseLit;
    falseLit.clear();

    bool True_confl = false;
    int i, j;
    CRef confl;

    for (i = 0, j = 0; i < c.size(); i++) {
        if (assignmentTrail.value(c[i]) == l_Undef) {
            propagationQueue.simpleEnqueue(~c[i]);
            c[j++] = c[i];
            confl = unitPropagator.simplePropagate();
            if (confl != CRef_Undef) break;
        } else {
            if (assignmentTrail.value(c[i]) == l_True){
                c[j++] = c[i];
                True_confl = true;
                confl = assignmentTrail.reason(var(c[i]));
                break;
            } else {
                falseLit.push(c[i]);
            }
        }
    }
    c.shrink(c.size() - j);

    if (confl != CRef_Undef || True_confl == true){
        simp_learnt_clause.clear();
        simp_reason_clause.clear();
        if (True_confl == true){
            simp_learnt_clause.push(c.last());
        }
        conflictAnalyzer.simpleAnalyze(confl, simp_learnt_clause, simp_reason_clause, True_confl, trailRecord);

        if (simp_learnt_clause.size() < c.size()){
            for (i = 0; i < simp_learnt_clause.size(); i++){
                c[i] = simp_learnt_clause[i];
            }
            c.shrink(c.size() - i);
        }
    }

    assignmentTrail.cancelUntilTrailSize(trailRecord);
}

template<int db_mark>
bool Solver::simplifyLearntDB() {
  int ci, cj;
    unsigned int nblevels;

    vec<CRef>& db = clauseDatabase.getDB<db_mark>();
    for (ci = 0, cj = 0; ci < db.size(); ci++) {
        CRef cr = db[ci];
        Clause& c = ca[cr];

        // Drop clause if already deleted
        if (ca[cr].mark() == 1) {
	  continue;
	}

        // Keep clause if already simplified
        if (c.simplified()) {
            db[cj++] = db[ci];
            continue;
        }

	
        // Attempt to simplify clause - drop clause if it is satisfied
        if (detachAndSimplify(cr)) continue;


        if (c.size() == 1) {
            // Enqueue and propagate unit clauses
            propagationQueue.enqueue(c[0]);
            if (unitPropagator.propagate() != CRef_Undef)
                return false;

            // delete the clause memory in logic
            c.mark(1);
            ca.free(cr);
        } else {
            // Add clause to appropriate clause database
            clauseDatabase.attachClause(cr);

            // Compute LBD of simplified clause
            nblevels = assignmentTrail.computeLBD(c);
            if (nblevels < static_cast<unsigned int>(c.lbd())) {
                c.set_lbd(nblevels);
            }

            // Try moving clause to core database
            if (db_mark != TIER2 || !clauseDatabase.upgradeToCore(cr, nblevels)) {
                db[cj++] = db[ci];
            }

            // Mark clause as simplified
            c.setSimplified(true);
        }
    }
    db.shrink(ci - cj);

    return true;
}

bool Solver::simplifyAll() {
    if (!ok || unitPropagator.propagate() != CRef_Undef)
        return ok = false;

    if (!simplifyLearntDB<CORE>()) return ok = false;

    if (!simplifyLearntDB<TIER2>()) return ok = false;
    // if (!simplifyLearntDB<LOCAL>()) return ok = false;

    clauseDatabase.checkGarbage();

    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// PROBLEM SPECIFICATION

bool Solver::addClause(vec<Lit>& ps) {
    assert(assignmentTrail.decisionLevel() == 0);

    // Only add clause if the solver is in a consistent state
    if (!ok) return false;

    // Make a copy of the input clause for proof output
    if (proofLogger.enabled()) ps.copyTo(add_oc);

    // Check if clause is satisfied and remove false/duplicate literals:
    sort(ps);
    Lit p; int i, j;
    for (i = j = 0, p = lit_Undef; i < ps.size(); i++) {
        if (assignmentTrail.value(ps[i]) == l_True || ps[i] == ~p)
            return true;
        else if (assignmentTrail.value(ps[i]) != l_False && ps[i] != p)
            ps[j++] = p = ps[i];
    }
    ps.shrink(i - j);

    // Output the simplified clause for proof output
    if (i != j) proofLogger.simplifyClause(add_oc, ps);

    // Mark solver as inconsistent if the clause is empty
    if (ps.size() == 0) return ok = false;
    
    // Add variable directly to trail if the clause is unit
    if (ps.size() == 1) {
        propagationQueue.enqueue(ps[0]);
        return ok = (unitPropagator.propagate() == CRef_Undef);
    }

    // Add the clause to the clause database
    clauseDatabase.addInputClause(ps);

    return true;
}

void Solver::writeClauseWithValues (const Clause& c, ostream& out){
  for (int i = 0; i < c.size(); ++i)
    out << c[i] << "(" << assignmentTrail.value(c[i]) << ", dl " << assignmentTrail.level(var(c[i])) << ") ";
}

void Solver::writeClauseWithValues (const vector<Lit>& c, ostream& out){
  for (uint i = 0; i < c.size(); ++i)
    out << c[i] << "(" << assignmentTrail.value(c[i]) << ", dl " << assignmentTrail.level(var(c[i])) << ") ";
}

void Solver::writeClauseWithValues (const vec<Lit>& c, ostream& out){
  for (int i = 0; i < c.size(); ++i)
    out << c[i] << "(" << assignmentTrail.value(c[i]) << ", dl " << assignmentTrail.level(var(c[i])) << ") ";
}
