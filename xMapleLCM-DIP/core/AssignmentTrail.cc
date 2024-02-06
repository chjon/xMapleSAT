/******************************************************************************[AssignmentTrail.cc]
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

#include "core/AssignmentTrail.h"
#include "core/Solver.h"

using namespace Minisat;

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

AssignmentTrail::AssignmentTrail(Solver& s)
    : ca(s.ca)
    , solver(s)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ACCESSORS

double AssignmentTrail::progressEstimate() const {
    double progress = 0;
    double F = 1.0 / nVars();

    for (int i = 0; i <= decisionLevel(); i++){
        int beg = i == 0 ? 0 : trail_lim[i - 1];
        int end = i == decisionLevel() ? nAssigns() : trail_lim[i];
        progress += pow(F, i) * (end - beg);
    }

    return progress / nVars();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// STATE MODIFICATION

void AssignmentTrail::assign(Lit p, CRef from) {
  genericAssign<false>(p, from);
}

void AssignmentTrail::simpleAssign(Lit p, CRef from){
  genericAssign<true>(p, from);
}

void AssignmentTrail::cancelUntilLevel(int level) {
    // Do nothing if the trail is already set at the correct level
    if (decisionLevel() <= level) return;

    // Backtrack
    cancelUntil<true>(trail_lim[level]);

    // Update current decision level
    trail_lim.shrink(trail_lim.size() - level);
}

void AssignmentTrail::cancelUntilTrailSize(int trailSize) {
    cancelUntil<false>(trailSize);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

template <bool simple>
inline void AssignmentTrail::genericAssign(Lit p, CRef from) {
    assert(value(p) == l_Undef);

    if (not simple) {
      // if (solver.conflicts == 2133) {
      // 	cout << "Set literal " << p << " due to ";
      // 	if (from == CRef_Undef) cout << " decision ";
      // 	else {
      // 		Clause& c = ca[from];
      // 		for (int i = 0; i < c.size(); ++i) cout << c[i] << "[" << value(c[i]) << ",dl " << level(var(c[i])) << "] ";
      // 	}
      // 	cout << " at DL " << decisionLevel() << endl;
      // }
#if DEBUG
      if (from != CRef_Undef) {
	Clause& c = ca[from];
	if (c.size() == 2) {
	  if (p == c[0]) {
	    assert(value(c[1]) == l_False);
	    assert(level(var(c[1])) == decisionLevel());
	  }
	  else {
	    assert(p == c[1]);
	    assert(value(c[0]) == l_False);
	    assert(level(var(c[0])) == decisionLevel());
	  }
	}
	else {
	  assert(c[0] == p);
	  bool allFalse = true;
	  int maxDL = -1;
	  for (int i = 1; i < c.size(); ++i) {
	    maxDL = max(maxDL,level(var(c[i])));
	    if (value(c[i]) != l_False) allFalse = false;
	  }
	  assert(allFalse);
	  assert(maxDL == decisionLevel());
	}
      }
#endif
    }
    // Notify event listeners
    if (!simple) solver.branchingHeuristicManager.handleEventLitAssigned(p, solver.conflicts);

    // Set variable value
    Var x = var(p);
    assigns[x] = lbool(!sign(p));
    if (simple) vardata[x].reason = from;
    else vardata[x] = VarData{from, decisionLevel()};
    trail.push_(p);
}

template<bool notifyListeners>
void AssignmentTrail::cancelUntil(int trailSize) {
    // Clear the values of the variables
    for (int c = trail.size() - 1; c >= trailSize; c--){
        Var x = var(trail[c]);
        assigns[x] = l_Undef;
	//	cout << "UnSet " << trail[c] << endl;

        // Notify event listeners
        if (notifyListeners) {
            solver.branchingHeuristicManager.handleEventLitUnassigned(
                trail[c],
                solver.conflicts,
                c > trail_lim.last()
            );
        }
    }

    // Decrease the size of the trail
    trail.shrink(trail.size() - trailSize);

    // Add remaining assignments to the queue
    solver.propagationQueue.batchEnqueue(trail, trailSize);
}
