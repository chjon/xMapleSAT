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
#include <set>
#include <ctime>
#include <fstream>
#include <iomanip>

using namespace Minisat;
using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////
// OPTIONS

static const char* _cat = "CORE";

static const char* _cat2 = "DIP";

static IntOption    opt_ccmin_mode             (_cat, "ccmin-mode",  "Controls conflict clause minimization (0=none, 1=basic, 2=deep)", 2, IntRange(0, 2));

static BoolOption   opt_learn_two_dip_clauses  (_cat2, "dip-2clauses",    "Learn two DIP clauses: UIP -> DIP and DIP -> conflict. If set to false, only DIP -> conflict is learned.", true);
static IntOption    opt_common_pair_DIP_min    (_cat2, "dip-pair-min",  "Specifies the minimum numer of times a DIP has to appear before we introduce it. (-1 means disabled)", 20, IntRange(-1, INT32_MAX));
static IntOption    opt_dip_type               (_cat2, "dip-type",  "Specifies the type of DIP computed (1 = middle, 2 = closest to conflict, 3 = random, 4 = heuristic)", 1, IntRange(1, 4));

static IntOption    opt_DIP_window_size         (_cat2, "dip-window-size",  "Introduce a DIP only if the sum of the activities of the pair is larger than the average of the last DIPs in a window of the given size (-1 means option disabled).", -1, IntRange(-1, INT32_MAX));

static IntOption    opt_DIP_max_LBD          (_cat2, "dip-glue",  "Introduce a DIP only if the clause LBD -> conflict has at most the given LBD  (-1 means option disabled).", -1, IntRange(-1, INT32_MAX));
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
  , dip_pair_threshold(opt_common_pair_DIP_min)
  , dip_window_size(opt_DIP_window_size)
  , dip_max_LBD(opt_DIP_max_LBD)
  , dip_type(opt_dip_type)
  , learn_two_DIP_clauses(opt_learn_two_dip_clauses)
    
    //////////////////////
    // Temporary variables
  , counter(0)

    /////////////
    // Statistics
  , max_literals(0)
  , tot_literals(0)
  , time_DIP(0)
  , conflicts_with_dip(0)
  , conflicts_with_dangerous_dip(0)
{
  if ( (dip_pair_threshold != -1 and
	dip_window_size != -1) or
       (dip_pair_threshold != -1 and
	dip_max_LBD != -1) or
       (dip_window_size != -1 and
	dip_max_LBD != -1)){
    cout << "ERROR: At most one of the options ,\"dip-pair-min\", \"dip-window\" and \"dip-glue\" can be enabled" << endl;
    exit(1);
  }  
}

/*_________________________________________________________________________________________________
|
|  analyze1UIP : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
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
void ConflictAnalyzer::analyze1UIP(CRef confl, vec<Lit>& out_learnt, int& out_btlevel, int& out_lbd) {
    // Generate conflict clause:
    mainLearnedClause.clear();
    getFirstUIPClause(confl, mainLearnedClause);
    max_literals += mainLearnedClause.size();

    // Simplify conflict clause:
    simplifyClause(out_learnt, mainLearnedClause);
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
    for (int j = 0; j < mainLearnedClause.size(); j++)
        seen[var(mainLearnedClause[j])] = false;
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
  |      * Returns whether a proper DIP has been found
  |  
  |________________________________________________________________________________________________@*/
bool ConflictAnalyzer::analyze (CRef confl, vec<Lit>& out_learnt, int& out_btlevel, int& out_lbd, vec<Lit>& out_learnt_UIP_to_DIP) {

  // First of all, depending on whether DIP learning has been performed, we compute
  // 1) [No-DIP learning]: mainLearnedClause is the 1UIP clause
  // 2) [DIP-learning]: mainLearnedClause is DIP -> conflict
  //                    out_learnt_UIP_to_DIP is 1UIP -> DIP
  mainLearnedClause.clear();
  bool found_DIP = getDIPLearntClauses(confl, mainLearnedClause, out_learnt_UIP_to_DIP);
  max_literals += mainLearnedClause.size();

  // Simplify conflict clause:
  simplifyClause(out_learnt, mainLearnedClause);
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
  // TODO (ALREADY FROM JONATHAN): can this be moved before updating the branching heuristic?
  // it is currently polluted by simplifyClause
  for (int j = 0; j < toClear.size(); j++)
    seen[toClear[j]] = false;
  toClear.clear();

  // Clear 'seen[]'
  for (int j = 0; j < mainLearnedClause.size(); j++) 
    seen[var(mainLearnedClause[j])] = false;

  assert(checkSeen());
  return found_DIP;
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
  // cout << endl << endl;
  // cout << "ENTER collectFirstUIP" << endl;
  assert(assignmentTrail.decisionLevel() > 0);
  int totalPathCount = 0;

  // cout << "Iterate over conflict: ";
  // for (int k = 0; k < ca[confl].size(); ++k) {
  //   Lit l = ca[confl][k];
  //   cout << l << " [" << assignmentTrail.value(l) << ",dl " << assignmentTrail.level(var(l)) << "] ";
  // }
  // cout << endl;

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

    // if (v == 14141) {
    //   cout << "We are at DL " << assignmentTrail.decisionLevel() << endl;
    //   cout << "assignmentTrail nVars " << assignmentTrail.nVars() << endl;
    //   cout << "reason of " << v << " " << assignmentTrail.reason(v) << endl;
    //   cout << "value of " << assignmentTrail.value(mkLit(v,true)) << endl;
    //   cout << "DL of " << v << " is " << assignmentTrail.level(v) << endl;
    //   cout << "Literal is " << p << endl;
    // }
    assert(assignmentTrail.level(v) > 0);
    Clause& rc = ca[assignmentTrail.reason(v)];
    enforceBinaryClauseInvariant(assignmentTrail, rc);

    // cout << "Iterate over reason for " << p << ": ";
    // for (int k = 0; k < rc.size(); ++k) {
    //   Lit l = rc[k];
    //   cout << l << "[" << assignmentTrail.value(l) << ", dl" << assignmentTrail.level(var(l)) << "] ";
    // }
    // cout << endl;

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

string graphVar2string(int v) {
  ++v;
  return to_string(v);    
}

string graphLit2String(Lit l) {
  Var v = var(l);
  bool isNegative = sign(l);
  string tmp = graphVar2string(v);
  if (isNegative) return "-" + tmp;
  else return tmp;
}

bool ConflictAnalyzer::ok_DIP (Lit dip1, Lit dip2, Lit UIP, CRef confl) {
  // Do a DFS starting from the "confl", only visiting notes at the current decision level, but not allowing to go
  // through dip1 or dip2. If we can reach UIP --> return false

  //cout << "DIP computation where dip1 " << dip1 << " dip2 " << dip2 << " UIP " << UIP << endl;
  
  vector<bool> visited(assignmentTrail.nVars(),false);
  visited[var(dip1)] = visited[var(dip2)] = true;

  assert(confl != CRef_Undef);
  stack<Lit> S;
  Clause& c = ca[confl];
  for (int i = 0; i < c.size(); ++i) {
    if (assignmentTrail.level(var(c[i])) == assignmentTrail.decisionLevel()) {
      S.push(c[i]);
      //cout << "Start with " << c[i] << " at DL " << assignmentTrail.level(var(c[i])) << endl;
    }    
  }

  while (not S.empty()) {
    Lit l = S.top(); S.pop();
    if (visited[var(l)]) continue;
    if (var(l) == var(UIP)) return false;
    visited[var(l)] = true;
    CRef r = assignmentTrail.reason(var(l));
    assert(r != CRef_Undef);
    Clause& reason = ca[r];
    for (int i = 0; i < reason.size(); ++i) {
      Lit l2 = reason[i];
      if (var(l2) == var(l)) continue;
      if (assignmentTrail.level(var(l2)) != assignmentTrail.decisionLevel()) continue;
      S.push(l2);
    }
  }
  return true;
}

bool ConflictAnalyzer::computeRandomDIP (int a, int b, TwoVertexBottlenecks& info, const DIPGraphEncoder& encoder, Lit& x, Lit& y) {
  // a == 1 if the first element of VertPairListA is an immediate ancestor of the conflict node
  // b == 1 if the first element of VertPairListB is an immediate ancestor of the conflict node
  
  vector<TwoVertexBottlenecks::VertPairInfo> listA = info.GetVertListA();
  int idxA;
  if (a == 1) {
    if (listA.size() == 1) return false;
    else idxA = rand()%(listA.size()-1) + 1;
  }
  else idxA = rand()%listA.size();
  
  x = encoder.Sam2Solver(listA[idxA].vertNum);

  vector<TwoVertexBottlenecks::VertPairInfo> listB = info.GetVertListB();  
  vector<Lit> candidatesY;
  for (int i = 0; i < listB.size(); ++i){
    if (b == 1 and i == 0) continue;
    if (listA[idxA].minPair <= listB[i].vertNum and
	listB[i].vertNum <= listA[idxA].maxPair)
      {
	if (listA[idxA].minAncestor <= listB[i].vertNum) { // skip ancestors (I do not think this is necessary)
	  //cout << "Skip candidate " << encoder.Sam2Solver(listB[i].vertNum) << " because it is an ancestor of " << x << endl;
	}
	else candidatesY.push_back(encoder.Sam2Solver(listB[i].vertNum));
      }
  }

  if (candidatesY.size() == 0) return false;
  y = candidatesY[rand()%candidatesY.size()];

  return true;
}

bool ConflictAnalyzer::computeClosestDIPToConflict (int a, int b, TwoVertexBottlenecks& info, const DIPGraphEncoder& encoder, Lit& x, Lit& y) {
  const vector<TwoVertexBottlenecks::VertPairInfo>& listA = info.GetVertListA();

  // a == 1 if the first element of VertPairListA is an immediate ancestor of the conflict node
  // b == 1 if the first element of VertPairListB is an immediate ancestor of the conflict node

  // Closest to the conflict
  int idxA;
  if (a == 1) {
    if (listA.size() == 1) return false;
    else idxA = 1;
  }
  else idxA = 0;

  const vector<TwoVertexBottlenecks::VertPairInfo>& listB = info.GetVertListB();    
  int idxB = -1;
  for (int i = 0; i < listB.size(); ++i) {
    if (b == 1 and i == 0) continue;
    if (listA[idxA].minPair <= listB[i].vertNum and
	listB[i].vertNum <= listA[idxA].maxPair and
	listA[idxA].minAncestor > listB[i].vertNum) { // avoid ancestors but I believe is not necessary
      idxB = i;
    }
  }

  if (idxB == -1) return false;

  x = encoder.Sam2Solver(listA[idxA].vertNum);
  y = encoder.Sam2Solver(listB[idxB].vertNum);

  return true;
}

bool ConflictAnalyzer::computeBestMiddleDIP (const TwoVertexBottlenecks& info, const DIPGraphEncoder& encoder, const vector<int>& pathA, const vector<int>& pathB, Lit& x, Lit& y) {

  if (pathA.size() <= 3 and pathB.size() <= 3) return false; // does not seem interesting
  
  
  const vector<TwoVertexBottlenecks::VertPairInfo>& listA = info.GetVertListA();
  const vector<TwoVertexBottlenecks::VertPairInfo>& listB = info.GetVertListB();

  // The idea is:
  // 1) Assume wlog that pathA is the longest path between pathA, pathB
  // 2) Among all DIPs in pathA, choose the one (DIP_A) that is closest to the middle position of pathA
  //    (i.e. ideally half the vertex in pathA)
  // 3) Having chosen DIP_A, choose DIP_B a vertex in pathB that, among the ones that form a DIP with DIP_A,
  //    is the one which is closes to the middle position inside pathB
  
  // Determine who is longest/shortest
  const vector<int>& longestPath = (pathA.size() > pathB.size() ? pathA : pathB);
  const vector<int>& shortestPath = (pathA.size() > pathB.size() ? pathB : pathA);
  const vector<TwoVertexBottlenecks::VertPairInfo>& longestList = (pathA.size() > pathB.size() ? listA : listB);
  const vector<TwoVertexBottlenecks::VertPairInfo>& shortestList = (pathA.size() > pathB.size() ? listB : listA);

  // First we compute the position in the path of each vertex of the lists (which are the DIP candidates)
  // For every k, posInLongestPath[k] = p  iff vertex longestList[k].vertNum appears in longestPath[p]
  // similarly for shortest
  vector<int> posInLongestPath(longestList.size());
  vector<int> posInShortestPath(shortestList.size());

  // Remember that pathA, pathB both start with the 1UIP (largest vertex) and finish with the conflict (vertex 0).
  // We know that both paths are strictly ordered from large to small
  // However, note that longest/shortest List are ordered strictly from smaller to larger
  for (int i = 0, j = longestPath.size() - 1; i < longestList.size(); ++i) {
    while (longestPath[j] != longestList[i].vertNum) --j;
    posInLongestPath[i] = j;
    //cout << "Pos of " << longestList[i].vertNum << " in longestPath is " << j << endl;
  }

  for (int i = 0, j = shortestPath.size() - 1; i < shortestList.size(); ++i) {
    while (shortestPath[j] != shortestList[i].vertNum) --j;
    posInShortestPath[i] = j;
    //cout << "Pos of " << shortestList[i].vertNum << " in shortestPath is " << j << endl;
  }

  int medianLongest = longestPath.size()/2;
  int medianShortest = shortestPath.size()/2;

  // We know that the longest path has size >= 3 (we have returned false otherwise)  
  int idxLongest = 0;
  int longestDistToMedian  = abs(medianLongest - posInLongestPath[0]);
  // Look for the one closest to the middle (DIP_A in the description)
  int k = 1;
  while (k < posInLongestPath.size() and longestDistToMedian > 0) {
    int dist = abs(medianLongest - posInLongestPath[k]);
    if (dist < longestDistToMedian) {
      longestDistToMedian = dist;
      idxLongest = k;
    }
    ++k;
  }

  //cout << "For longest we choose DIP " << longestList[idxLongest].vertNum << " at position " << posInLongestPath[idxLongest] << endl;

  // minV and maxV tell that any number inside [minV, maxV] is DIP-compatible with DIP_A
  int minV = longestList[idxLongest].minPair;
  int maxV = longestList[idxLongest].maxPair;

  int idxShortest = -1;
  int shortestDistToMedian = 1e9;
  
  k = 0;
  while (k < posInShortestPath.size() and shortestDistToMedian > 0) {
    if (shortestList[k].vertNum < minV or shortestList[k].vertNum > maxV) {++k; continue;}
    int dist = abs(medianShortest - posInShortestPath[k]);
    if (dist < shortestDistToMedian) {
      shortestDistToMedian = dist;
      idxShortest = k;
    }
    ++k;
  }

  assert(idxShortest >= 0);
  // if (idxShortest < 0) {  // Just in case. If assertion fails sometime, add this code
  //   return false; 
  // }

  //cout << "For shortest we choose DIP " << shortestList[idxShortest].vertNum << " at position " << posInShortestPath[idxShortest] << endl;


  if (posInShortestPath[idxShortest] == shortestPath.size() - 2 and
      posInLongestPath[idxLongest] == longestPath.size() - 2)
    {
      //cout << "TWO IMMEDIATE ANCESTORS" << endl;
      return false;}
  
  
  x = encoder.Sam2Solver(longestList[idxLongest].vertNum);
  y = encoder.Sam2Solver(shortestList[idxShortest].vertNum);


  return true;

}

// Returns whether a DIP has been found
// INPUT: a, b (return values of DIP computation)
//        info --> class that computed all DIPs
//        encoder --> mapping between ints in the DIP computation algorithm and literals
//        pathA, pathB --> paths in the DIP computation
// OUTPUT: x, y --> DIP pair
//
// I would recommend to have a look at computeBestMiddleDIP. It uses all these data structures for computing
// another type of DIP
bool ConflictAnalyzer::computeHeuristicDIP (int a, int b, const TwoVertexBottlenecks& info, const DIPGraphEncoder& encoder, const vector<int>& pathA, const vector<int>& pathB, Lit& x, Lit& y) {
  const vec<double>& acts = branchingHeuristicManager.getActivity(); // indexed by variable
  
  
  return true;
}

int ConflictAnalyzer::computeLBD_DIP2Conflict (CRef confl, Lit x, Lit y) {

  int dipReached = 0;
  Lit p = lit_Undef;
  int index = assignmentTrail.nAssigns() - 1;
  int heightX = -1, heightY = -1;
  vector<Lit> afterLits;
  do {
    assert(confl != CRef_Undef); // (otherwise should be UIP)
    Clause& c = ca[confl];

    // For binary clauses, we don't rearrange literals in propagate(), so check and make sure the first is an implied lit.
    if (p != lit_Undef && c.size() == 2 && assignmentTrail.value(c[0]) == l_False) {
      assert(assignmentTrail.value(c[1]) == l_True);
      std::swap(c[0], c[1]);
    }
	
    // Iterate through every literal that participates in the conflict graph
    for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++) {
      Lit q = c[j];
      if (seen3[var(q)] || assignmentTrail.level(var(q)) == 0) continue;

      // Mark variable as seen      
      seen3[var(q)] = 1;

      // Increment the number of paths if the variable is assigned at the decision level that resulted in conflict
      if (assignmentTrail.level(var(q)) < assignmentTrail.decisionLevel())
	afterLits.push_back(q);
    }
        
    // Select next clause to look at:
    while (dipReached != 2) {
      while (!seen3[var(assignmentTrail[index--])]);
      p = assignmentTrail[index+1];
      if (p == x) heightX = index+1;
      else if (p == y) heightY = index+1;	  
      if (p != x and p != y) break;
      else dipReached++;
    }
	
    p     = assignmentTrail[index+1];
    confl = assignmentTrail.reason(var(p));
	
    // Mark variable as unseen: it is either at or after the first UIP
    seen3[var(p)] = 0;
    
  } while (dipReached != 2);

  static vector<int> levels(assignmentTrail.nVars()+1,0);
  static int times = 0;
  ++times;
  if (times > 1e9) {
    for (int& z : levels) z = 0;
    times = 1;
  }

  int LBD = 0;
  for (auto z : afterLits) {
    int lev = assignmentTrail.level(var(x));
    if (levels[lev] != times) {levels[lev] = times; ++LBD;}    
  }
  
  // This was the idea of forcing the heuristic to immediately branch on one of the DIPs
  //solver.branchingHeuristicManager.nextDecisionBestOf(~x,~y); 

  // Clear marks
  for (int i = 0; i < afterLits.size(); ++i) seen3[var(afterLits[i])] = false;
  seen3[var(x)] = seen3[var(y)] = false;

  assert(checkSeen3());
  return LBD + 1;
}

bool ConflictAnalyzer::computeDIPClauses (int a, int b, CRef confl, TwoVertexBottlenecks& info, const DIPGraphEncoder& encoder, vec<Lit>& clause_to_learn, vec<Lit>& clause_to_learn2, Lit UIP, vector<int>& pathA, vector<int>& pathB) {
  assert(clause_to_learn.size() == 0);

  CRef originalConflict = confl;
    
  Lit x, y; // {x,y} are a DIP
  // 3 DIP computations: random, closest to conflict, the one in the middle
  
  bool dip_found = false;
  if (dip_type == MIDDLE_DIP) dip_found = computeBestMiddleDIP(info,encoder,pathA,pathB,x,y);
  else if (dip_type == CLOSEST_TO_CONFLICT) dip_found = computeClosestDIPToConflict(a,b,info,encoder,x,y);
  else if (dip_type == RANDOM_DIP) dip_found = computeRandomDIP(a,b,info,encoder,x,y);
  else if (dip_type == HEURISTIC_DIP) dip_found = computeHeuristicDIP(a,b,info,encoder,pathA,pathB,x,y);
  else assert(false);

  if (not dip_found) return false;
  
  // //if (not computeRandomDIP(a,b,info,encoder,x,y)) return false;
  // //if (not computeClosestDIPToConflict(a,b,info,encoder,x,y)) return false;
  // if (not computeBestMiddleDIP(info,encoder,pathA,pathB,x,y)) return false;
  
  //cout << "DIP " << x << " " << y << " at conflict " << solver.conflicts << endl;


  if (dip_pair_threshold != -1) { // common pair used
    solver.branchingHeuristicManager.notifyDIPCandidateCommonPair(~x,~y); 
    if (not solver.branchingHeuristicManager.isPairCommon(~x,~y,dip_pair_threshold)) return false;
  }
  else if (dip_window_size != -1) { // use sliding window
    solver.branchingHeuristicManager.notifyDIPCandidateWindow(~x,~y);
    if (not solver.branchingHeuristicManager.isPairBetterThanAverage(~x,~y)) return false;
  }
  // else if (dip_max_LBD != -1) { // use maxLBD for DIP clause
  //   int expectedLBD = computeLBD_DIP2Conflict(confl,x,y) + 1;
  //   if (expectedLBD > dip_max_LBD) return false;
  // }


  // 1) Compute first clause: DIP ^ after -> conflict
  // This is done starting from conflict, doing standard 1UIP reasoning but stopping as soon as we hit the
  // two elements of the DIP
  assert(checkSeen3());
  vector<Lit> toIncreaseActivity;
  int dipReached = 0;
  Lit p = lit_Undef;
  int index = assignmentTrail.nAssigns() - 1;
  int heightX = -1, heightY = -1;
  vector<Lit> afterLits;
  do {
    assert(confl != CRef_Undef); // (otherwise should be UIP)
    Clause& c = ca[confl];

    // For binary clauses, we don't rearrange literals in propagate(), so check and make sure the first is an implied lit.
    if (p != lit_Undef && c.size() == 2 && assignmentTrail.value(c[0]) == l_False) {
      assert(assignmentTrail.value(c[1]) == l_True);
      std::swap(c[0], c[1]);
    }
	
    // Iterate through every literal that participates in the conflict graph
    for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++) {
      Lit q = c[j];
      if (seen3[var(q)] || assignmentTrail.level(var(q)) == 0) continue;

      // Increment activity heuristic
      //branchingHeuristicManager.handleEventLitInConflictGraph(q, solver.conflicts);
      toIncreaseActivity.push_back(q);
      
      // Mark variable as seen      
      seen3[var(q)] = 1;

      // Increment the number of paths if the variable is assigned at the decision level that resulted in conflict
      if (assignmentTrail.level(var(q)) < assignmentTrail.decisionLevel())
	afterLits.push_back(q);
    }
        
    // Select next clause to look at:
    while (dipReached != 2) {
      while (!seen3[var(assignmentTrail[index--])]);
      p = assignmentTrail[index+1];
      if (p == x) heightX = index+1;
      else if (p == y) heightY = index+1;	  
      if (p != x and p != y) break;
      else dipReached++;
    }
	
    p     = assignmentTrail[index+1];
    confl = assignmentTrail.reason(var(p));
	
    // Mark variable as unseen: it is either at or after the first UIP
    seen3[var(p)] = 0;
    
  } while (dipReached != 2);
  
  // This was the idea of forcing the heuristic to immediately branch on one of the DIPs
  //solver.branchingHeuristicManager.nextDecisionBestOf(~x,~y); 

  // Clear marks
  for (int i = 0; i < afterLits.size(); ++i) seen3[var(afterLits[i])] = false;
  seen3[var(x)] = seen3[var(y)] = false;
  assert(checkSeen3());

  if (dip_max_LBD != -1) {
    int LBD = assignmentTrail.computeLBD(afterLits) + 1;
    if (LBD > dip_max_LBD) return false;
  }
  
  for (Lit z : toIncreaseActivity)
    branchingHeuristicManager.handleEventLitInConflictGraph(z, solver.conflicts);
  
  // 2) add definition clauses
  auto [newDef, extLit] = erManager->addDIPExtensionVariable(~x,~y);

  // cout << extLit << " <-> " << ~x << " v "  << ~y << endl;
  // cout << ~extLit << " <-> " << x << " ^ "  << y << endl;
  if (not newDef) {
    // This is for debugging
    // cout << "Careful: this is not a new definition" << endl;
    // cout << x << " --> " << assignmentTrail.value(x) << "[dl" << assignmentTrail.level(var(x)) << "]" << endl;
    // cout << y << " --> " << assignmentTrail.value(y) << "[dl" << assignmentTrail.level(var(y)) << "]" << endl;
    // cout << extLit << " --> " << assignmentTrail.value(extLit) << "[dl" << assignmentTrail.level(var(extLit)) << "]" << endl;

    // Here we compute the truth value of the conjuntion x ^ y (right) and also of the negation of the
    // extension variable (left). If they do not match we do not use DIP
    // (TO BE FIXED, since this happens quite often)
    
    lbool right = l_Undef;
    if (assignmentTrail.value(x) == l_True and assignmentTrail.value(y) == l_True) right = l_True;
    else if (assignmentTrail.value(x) == l_False or assignmentTrail.value(y) == l_False) right = l_False;
    lbool left = assignmentTrail.value(~extLit);
    if ((right == l_True and left == l_False) or
	(right == l_False and left == l_True)) {
      ++conflicts_with_dangerous_dip;
      // static int times = 0;
      // ++times;
      // if (times % 1000 == 1) {
      // 	cout << "CAREFUL: the definition is falsified in conflict " << solver.conflicts << " at DL " << assignmentTrail.decisionLevel() << " (has happened " << times << " times)" << endl;
      // }
      //      exit(1);
      return false;
    }
  }

  if (newDef and solver.produce_proof) { // Add definition to proof
    static vec<Lit> ERClause;
    ERClause.clear();
    ERClause.push(~extLit);     ERClause.push(~x);     ERClause.push(~y); 
    solver.proofLogger.addClause(ERClause);

    ERClause.clear();
    ERClause.push(extLit);    ERClause.push(x);     
    solver.proofLogger.addClause(ERClause);
    
    ERClause.clear();
    ERClause.push(extLit);    ERClause.push(y);     
    solver.proofLogger.addClause(ERClause);
  }

  clause_to_learn.push(extLit);
  for (auto x : afterLits) clause_to_learn.push(x);

  if (not learn_two_DIP_clauses) return true;
    


  // 3) UIP ^ before --> DIP
  // This is done as standard 1UIP but starting with the two DIP lits being marked
  // Also, we start having a look at the trail at the max height of the DIPs lits  

  assert(assignmentTrail[heightX] == x and assignmentTrail[heightY] == y);

  int open = 1; // not two because in fact we already take the reason of the highest DIP lit. Hence
                // this lit has already been "regressed"
  p = (heightX > heightY ? x : y);
  index = max(heightX,heightY) - 1;
  vector<Lit> beforeLits;
  confl = assignmentTrail.reason(var(p));
  seen3[var(x)] = seen3[var(y)] = true;
  // cout << "Començo amb x " << x << " alçada " << heightX << endl;
  // cout << "L'altre es y " << y << " alçada " << heightY << endl;
  do {
    assert(confl != CRef_Undef); // (otherwise should be UIP)
    Clause& c = ca[confl];

    // For binary clauses, we don't rearrange literals in propagate(), so check and make sure the first is an implied lit.
    if (p != lit_Undef && c.size() == 2 && assignmentTrail.value(c[0]) == l_False) {
      assert(assignmentTrail.value(c[1]) == l_True);
      std::swap(c[0], c[1]);
    }

    // Iterate through every literal that participates in the conflict graph
    for (int j = 1; j < c.size(); j++) {
      Lit q = c[j];
      if (seen3[var(q)] || assignmentTrail.level(var(q)) == 0) continue;

      // Increment activity heuristic
      branchingHeuristicManager.handleEventLitInConflictGraph(q, solver.conflicts);
      // Mark variable as seen      
      seen3[var(q)] = 1;
      
      // Increment the number of paths if the variable is assigned at the decision level that resulted in conflict
      if (assignmentTrail.level(var(q)) < assignmentTrail.decisionLevel())
	beforeLits.push_back(q);
      else ++open;
    }
        
    // Select next clause to look at:
    while (!seen3[var(assignmentTrail[index--])]);
    p = assignmentTrail[index+1];
    confl = assignmentTrail.reason(var(p));
    //	cout << "Ara trec la reason de " << p << endl;
	
    // Mark variable as unseen: it is either at or after the first UIP
    seen3[var(p)] = 0;
    --open;
	
  } while (open > 0);


  clause_to_learn2.push(~p);
  clause_to_learn2.push(~extLit); // Two first lits are of highest DL
  // First is UIP
  assert(assignmentTrail.level(var(clause_to_learn2[0])) == assignmentTrail.decisionLevel());  
  for (auto l : beforeLits) {
    clause_to_learn2.push(l);
    assert(assignmentTrail.level(var(l)) < assignmentTrail.decisionLevel());
  }

  // Clear marks
  for (int i = 0; i < clause_to_learn2.size(); ++i) seen3[var(clause_to_learn2[i])] = false;
  seen3[var(p)] = seen3[var(x)] = seen3[var(y)] = false;
  
  assert(checkSeen3());
  
  return true;
}

void writeFinalInfoGraph (ofstream& out, set<Lit>& current, set<Lit>& previous, Lit UIP, bool writePrevious) {
  for (auto lit : current)
    if (lit != UIP) 
      out << "\"" << graphLit2String(lit) << "\"[ fillcolor = lightblue, style = filled ];" << endl;  
  
  if (writePrevious)
    for (auto lit : previous)
      out << "\"" << graphLit2String(lit) << "\"[ fillcolor = lightgrey, style = filled ];" << endl;
  
  out << "\"" << graphLit2String(UIP) << "\"[ fillcolor = red, style = filled ];" << endl;
  out << "}" << endl;
}

void ConflictAnalyzer::writeEdgeInGraph (ofstream& out, const Lit& orig, const Lit& dest, bool colored) {
  if (var(dest) == assignmentTrail.nVars()) out << "\"" << graphLit2String(orig) << "\" -> CONFLICT";
  else out << "\"" << graphLit2String(orig) << "\" -> \"" << graphLit2String(dest) << "\" ";
  if (colored)  out << "[color=red];";
  out << endl;
}

void ConflictAnalyzer::writeDIPComputationInfo (TwoVertexBottlenecks& info, DIPGraphEncoder& encoder, const vector<int>& predecessors, const vector<Lit>& predecessorsLits, const vector<int>& predIndex, const vector<Lit>& literalsInAnalysis, bool foundDIP, const vector<int>& pathA, const vector<int>& pathB) {
  cout << endl;
  cout << endl;
  cout << string(60,'=') << endl;
  cout << "Conflict " << solver.conflicts << " at level " << assignmentTrail.decisionLevel() << endl;
  cout << "N " << encoder.numVertices() << endl;
  cout << "predecessors: ";
  for (int x : predecessors) cout << x << " ";
  cout << endl;
  cout << "predIndex: ";
  for (int x : predIndex) cout << x << " ";
  cout << endl;
  
  cout << "Literals in analysis: ";
  for (auto lit : literalsInAnalysis) cout << lit << " ";
  cout << endl;
  cout << "Encoding: " << endl;
  for (auto lit : predecessorsLits)
    cout << "Lit " << lit << " is DIP algorithm vertex " << encoder.Solver2Sam(lit) << endl;
  
  cout << "Found DIP: " << (foundDIP?"true":"false") << endl;
  
  vector<TwoVertexBottlenecks::VertPairInfo> listA = info.GetVertListA();
  vector<TwoVertexBottlenecks::VertPairInfo> listB = info.GetVertListB();

  cout << string(30,'-') << endl;
  cout << "pathA: ";
  for (int x : pathA) cout << x << " ";
  cout << endl;
  cout << "pathB: ";
  for (int x : pathB) cout << x << " ";
  cout << endl << endl;

  cout << "List A: " << endl;
  for (int i = 0; i < listA.size(); ++i){
    cout << "Vertex[" << i << "]: " << listA[i].vertNum << " --> lit " << encoder.Sam2Solver(listA[i].vertNum) << endl;
    cout << "minPair: " << listA[i].minPair << endl; 
    cout << "maxPair: " << listA[i].maxPair << endl;
    cout << "minAncestor: " << listA[i].minAncestor << endl;
  }

  cout << endl;
  cout << "List B: " << endl;
  for (int i = 0; i < listB.size(); ++i){
    cout << "Vertex[" << i << "]: " << listB[i].vertNum << " --> lit " << encoder.Sam2Solver(listB[i].vertNum) << endl;
    cout << "minPair: " << listB[i].minPair << endl; 
    cout << "maxPair: " << listB[i].maxPair << endl;
    cout << "minAncestor: " << listB[i].minAncestor << endl;
  }
}
						
inline bool ConflictAnalyzer::getDIPLearntClauses (CRef confl, vec<Lit>& out_learnt, vec<Lit>& out_learnt_UIP_to_DIP) {
  // Stuff for writing the conflict graph (useful for debugging/understanding)
  // Currently the writing is interleaved with the 1UIP computation
  // It might be cleaner to do it separately (but I think that is not a problem for efficiency)
  string filename1 = "graph-"+to_string(solver.conflicts) + ".dot"; //conflict graph
  string filename2 = "graph-current-"+to_string(solver.conflicts) + ".dot"; // conflict graph with only last-DL lits
  ofstream outAll;
  ofstream outCurrent;
  set<Lit> currentNodes;
  set<Lit> previousNodes;

  //if (solver.conflicts == 45) exit(1);
  //bool write = solver.conflicts == 44; // To write conflict graph only of a concrete conflict
#define write 0 // quicker for release mode
  if (write) {
    outAll.open(filename1.c_str(),fstream::out);
    outCurrent.open(filename2.c_str(),fstream::out);
    outAll << "digraph D {" << endl;
    outCurrent << "digraph D {" << endl;
  }

  // Data structure for generating correct input for DIP-detection algorithm
  DIPGraphEncoder encoder(assignmentTrail.nVars());
  vector<Lit> predecessorsLits; // only needed for DIP-detection algorithm
  vector<int> predIndex; // only needed for DIP-detection algorithm
  vector<Lit> literalsInAnalysis; // only needed for DIP-detection algorithm
  // The latter stores the literals as they have appeared in the analysis (inverse topological order)
  // This will allow us to assign numbers respecting the order (the algorithm needs this)
  // Remember that conflict node is the sink (node 0), and the 1UIP is the source (node with largest N)
  // Once the analysis is done, we only have to convert the predecessorLits into a vector<int>
  // to get the correct input to DIP-detection


  CRef origConfl = confl; // needed because confl is modified during the 1UIP detection    

  vector<Lit> litsToNotifyBranchingHeuristic; // ALBERT: should explain and understand this

  // Start by adding the edges from the negation of the lits  the conflicting clause to the "CONFLICT" node
  literalsInAnalysis.push_back(mkLit(assignmentTrail.nVars(),false)); // Fake literal corresponding to conflict
  predIndex.push_back(predecessorsLits.size());      
  Clause& confClause = ca[confl];  
  for (int i = 0; i < confClause.size(); ++i) {
    if (write) writeEdgeInGraph(outAll,~confClause[i],mkLit(assignmentTrail.nVars(),false),confClause.learnt());
    if (assignmentTrail.level(var(confClause[i])) == assignmentTrail.decisionLevel()) {
      predecessorsLits.push_back(~confClause[i]);       	
      if (write) {
	currentNodes.insert(~confClause[i]);
	writeEdgeInGraph(outCurrent,~confClause[i],mkLit(assignmentTrail.nVars(),false),confClause.learnt());
      }
    }
    else if (write) previousNodes.insert(~confClause[i]);
  }

  // Start standard 1UIP learning algorithm
  // We generate the 1UIP clause because if there is no DIP, we are going to learn 1UIP
  Lit UIP = lit_Undef;
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
      //branchingHeuristicManager.handleEventLitInConflictGraph(q, solver.conflicts);
      litsToNotifyBranchingHeuristic.push_back(q);
      // We do not increase the activity score here to give some flexibility to use other approaches (e.g.
      // only increase the ones from DIP to conflict?)
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
    if (pathC > 0) {
      Clause& cTmp = ca[confl];
      literalsInAnalysis.push_back(p);
      predIndex.push_back(predecessorsLits.size());	    
      for (int i = 0; i < cTmp.size(); ++i)
	if (var(cTmp[i]) != var(p)) {
	  if (write) writeEdgeInGraph(outAll,~cTmp[i],p,cTmp.learnt());		
	  if (assignmentTrail.level(var(cTmp[i])) == assignmentTrail.decisionLevel()) {
	    predecessorsLits.push_back(~cTmp[i]);				
	    if (write) {
	      currentNodes.insert(~cTmp[i]);
	      writeEdgeInGraph(outCurrent,~cTmp[i],p,cTmp.learnt());
	    }
	  }
	  else if (write) previousNodes.insert(~cTmp[i]);
	}
    }
    else if (pathC == 0) { // 1UIP found
      predIndex.push_back(predecessorsLits.size());
      literalsInAnalysis.push_back(p);
      UIP = p;
    }

  } while (pathC > 0);

  if (write) {
    writeFinalInfoGraph(outAll,currentNodes,previousNodes,UIP,true);
    writeFinalInfoGraph(outCurrent,currentNodes,previousNodes,UIP,false);
    outAll.close();
    outCurrent.close();
  }

  bool foundDIP = false;

  clock_t start = clock(), end;
  
  vector<int> predecessors;
  for (auto lit : literalsInAnalysis) encoder.Solver2Sam(lit); // encode in topological order
  for (auto lit : predecessorsLits) predecessors.push_back(encoder.Solver2Sam(lit)); // map vec<Lit> to vec<int>
  TwoVertexBottlenecks dip;
  vector<int> pathA, pathB; // the two disjoint paths
  int res = dip.CalcBottlenecks(predecessors,predIndex,pathA,pathB);
  
  foundDIP =  (res > 0);
  
  
  if (write) 
    writeDIPComputationInfo(dip,encoder,predecessors,predecessorsLits,predIndex,literalsInAnalysis,foundDIP,pathA,pathB);
  
  if (foundDIP) {
    
    ++conflicts_with_dip;
    
    // Decode
    //    The return code res will equal:
    //                  res = 4 + 2a + b
    //    where:  a==1 if the first element of VertPairListA is an immediate
    //                 ancestor of the sink node
    //            b==1 if the first element of VertPairListB is an immediate
    //                 ancestor of the sink node
    //    and a, b are otherwise equal to 0.    
    res -=4;
    int b = res%2;
    res/=2;
    int a = res;     
    
    vec<Lit> dip_clause_to_learn;
    vec<Lit> dip_clause_to_learn2;
    
    bool ok = computeDIPClauses(a,b,origConfl,dip,encoder,dip_clause_to_learn,dip_clause_to_learn2,UIP,pathA, pathB);
    // If DIP-clause learning seems a good idea
    if (ok) {
	// Clear seen marks
      for (int i = 1; i < out_learnt.size(); ++i) seen[var(out_learnt[i])] = false; // clearing the 1st UIP lemma
      if (learn_two_DIP_clauses) for (int i = 0; i < dip_clause_to_learn.size(); ++i) seen[var(dip_clause_to_learn[i])] = true;
#if DEBUG
      for (int k = 1; k < dip_clause_to_learn.size(); ++k)
	assert(assignmentTrail.value(dip_clause_to_learn[k]) == l_False);
#endif 	  
      // Copy dip_clause_to_learn (DIP -> conflict) to out_learnt
      out_learnt.clear();
      for (int i = 0; i < dip_clause_to_learn.size() ;++i)
	out_learnt.push(dip_clause_to_learn[i]);
      
      if (learn_two_DIP_clauses) {
	// Copy dip_clause_to_learn2 (1UIP -> conflict) to out_learnt_UIP_to_DIP
	out_learnt_UIP_to_DIP.clear();
	for (int i = 0; i < dip_clause_to_learn2.size() ;++i)
	  out_learnt_UIP_to_DIP.push(dip_clause_to_learn2[i]);
      }
      
      end = clock();
      time_DIP += double(end - start)/CLOCKS_PER_SEC;
      return true;
    }
    else { // DIP existed but not good enough or some generated clause is dangerous (to be fixed?)
      // In this case learn 1UIP clause (which is inside out_learnt)
      end = clock();
      time_DIP += double(end - start)/CLOCKS_PER_SEC;
      
      out_learnt[0] = ~p; 
      seen[var(p)] = true;
      for (auto l : litsToNotifyBranchingHeuristic) // If 1UIP learning, activity should be like 1UIP learning
	branchingHeuristicManager.handleEventLitInConflictGraph(l, solver.conflicts);
      return false; // no DIP found	
    }
  }
  end = clock();    
  time_DIP += double(end - start)/CLOCKS_PER_SEC;
  
  
  //    out << "UIP is " << p << endl;
  // Add first UIP literal at index 0
  for (auto l : litsToNotifyBranchingHeuristic)
    branchingHeuristicManager.handleEventLitInConflictGraph(l, solver.conflicts);

  out_learnt[0] = ~p;
  seen[var(p)] = true;

  return false; // no DIP found
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

void ConflictAnalyzer::notifyERManager(ERManager* erm) {
  erManager = erm;
}

bool ConflictAnalyzer::checkSeen ( ){
  int t = 0;
  for (int i = 0; i < seen.size(); ++i)
    if (seen[i] and assignmentTrail.value(mkLit(i)) != l_Undef and assignmentTrail.level(i) != 0) {
      cout << "Error with " << i << endl;
      ++t;
    }
  return t == 0;
}

bool ConflictAnalyzer::checkSeen3 ( ){
  int t = 0;
  for (int i = 0; i < seen3.size(); ++i)
    if (seen3[i] and assignmentTrail.value(mkLit(i)) != l_Undef and assignmentTrail.level(i) != 0) {
      cout << "Error3 with " << i << endl;
      ++t;
    }
  return t == 0;
}
