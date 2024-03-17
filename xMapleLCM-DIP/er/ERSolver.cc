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
// OPTIONS
static const char* _cat2 = "DIP";
static BoolOption   opt_compute_dip            (_cat2, "compute-dip",   "Compute DIP.", true);
static BoolOption   disabling_dip              (_cat2, "disabling-dip",   "Allows dynamically disabling DIP computation if the perc. of decisions on exteded variables is low.", false);

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

ERSolver::ERSolver()
  : Solver()
  , erManager(*this)
  , use_dip(opt_compute_dip)
  , allow_dip_disabling(disabling_dip)
{
  this->conflictAnalyzer.notifyERManager(&erManager);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

void ERSolver::disable_DIP_computation_if_needed ( ){
  if (not allow_dip_disabling) return;
  static bool some_very_high = false;
  static bool some_evaluation = false;
  static vector<double> window;  
  if (use_dip and not some_very_high) {
    if (conflicts%10000 == 0) {
      window.push_back(double(erManager.branchOnExt)/branchingHeuristicManager.decisions*100);
      if ((not some_evaluation and window.size() == 20) or (some_evaluation and window.size() == 10)) {
	some_evaluation = true;
	//	cout << "Window full" << endl;
	if (window.back() >= 5) {/*cout << "Some very high" << endl;*/ some_very_high = true;}
	if (some_very_high or
	    (window.back() - window[0] > 0)) {/*cout << "Continue" << endl;*/}
	else {/*cout << "Deactivate" << endl;*/use_dip = false;}
	window.clear();
      }
    }
  }
}

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
  // ALBERT erManager.checkGenerateDefinitions(conflicts);
#endif 

#if ER_USER_ADD_LOCATION == ER_ADD_LOCATION_AFTER_RESTART
  // Add extension variables
  // ALBERT erManager.introduceExtVars();
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

      disable_DIP_computation_if_needed();
      // if (conflicts % 10000 == 0)
      // 	cout << "We had " << branchingHeuristicManager.decisions << " decisions of which only " << double(erManager.branchOnExt)/branchingHeuristicManager.decisions*100 << " were on extended" << endl;

      
      // CONFLICT
      assert(conflictAnalyzer.checkSeen());
      conflicts++;

      if (assignmentTrail.decisionLevel() == 0) return l_False;

      clauseDatabase.handleEventConflicted(conflicts);
      branchingHeuristicManager.handleEventConflicted(confl, conflicts);
      learnt_clause.clear();
      learnt_clause_UIP_to_DIP.clear();

      if (use_dip) {     
	if (conflictAnalyzer.analyze(confl, learnt_clause, backtrack_level, lbd, learnt_clause_UIP_to_DIP))
	  dip_conflicts++;
      }
      else conflictAnalyzer.analyze1UIP(confl, learnt_clause, backtrack_level, lbd);
	
      assert(conflictAnalyzer.checkSeen());
      assignmentTrail.cancelUntilLevel(backtrack_level);

      if (use_dip) {
	// EXTENDED RESOLUTION - substitute disjunctions with extension variables. This must be
	// called after backtracking because extension variables might need to be propagated.
	bool subst = erManager.substitute(learnt_clause);
	
	
	if (subst) {
	  // Check for largest DL not counting UIP and move it to 2n pos
	  int max_pos = 1;
	  int max_level = assignmentTrail.level(var(learnt_clause[1]));
	  for (int k = 2; k < learnt_clause.size(); ++k)
	    if (assignmentTrail.level(var(learnt_clause[k])) > max_level) {
	      max_level = assignmentTrail.level(var(learnt_clause[k]));
	      max_pos = k;
	    }
	  
	  if (max_pos != 1)
	    swap(learnt_clause[1],learnt_clause[max_pos]);
	  
	  assignmentTrail.cancelUntilLevel(max_level);
	  lbd = assignmentTrail.computeLBD(learnt_clause); // it might have changed
	}
      }
	    
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
      // ALBERT erManager.checkGenerateDefinitions(conflicts); DISABLED FOR DIP-learning
#endif

#if ER_USER_ADD_LOCATION == ER_ADD_LOCATION_AFTER_CONFLICT
      // Add extension variables
      // ALBERT erManager.introduceExtVars(); DISABLED for DIP-learning
#endif

#if ER_ENABLE_GLUCOSER
      cout << "Glucoser" << endl; exit(1); // Should never enter here in DIP-learning
      erManager.generateLER();
      erManager.introduceExtVars(ERManager::HeuristicType::LER);
#endif

      // Output to proof file
      proofLogger.addClause(learnt_clause);	    
	    
      if (learnt_clause_UIP_to_DIP.size()) {
	assert(learn_two_DIP_clauses);
	// It is not obvious whether we should decrement the LBD by 1 or not!
	cr == clauseDatabase.addLearntClause(learnt_clause_UIP_to_DIP, assignmentTrail.computeLBD(learnt_clause_UIP_to_DIP) - 1, conflicts);

	// Should we do that????
	//restartHeuristicManager.handleEventLearntClause(lbd);

	proofLogger.addClause(learnt_clause_UIP_to_DIP);
      }
	    	    

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



