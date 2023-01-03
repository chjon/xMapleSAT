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
    : variableDatabase(s.variableDatabase)
    , ca(s.ca)
    , solver(s)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ACCESSORS

double AssignmentTrail::progressEstimate() const {
    double progress = 0;
    double F = 1.0 / variableDatabase.nVars();

    for (int i = 0; i <= decisionLevel(); i++){
        int beg = i == 0 ? 0 : trail_lim[i - 1];
        int end = i == decisionLevel() ? nAssigns() : trail_lim[i];
        progress += pow(F, i) * (end - beg);
    }

    return progress / variableDatabase.nVars();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// STATE MODIFICATION

void AssignmentTrail::assign(Lit p, CRef from) {
    assert(variableDatabase.value(p) == l_Undef);
    const Var v = var(p);
    solver.branchingHeuristicManager.handleEventLitAssigned(p, solver.conflicts);
    variableDatabase.setVar(v, lbool(!sign(p)));
    vardata[v] = VarData{from, decisionLevel()};
    trail.push_(p);
}

void AssignmentTrail::cancelUntil(int level) {
    // Do nothing if the trail is already set at the correct level
    if (decisionLevel() <= level) return;

    // Clear the values of the variables
    for (int c = trail.size() - 1; c >= trail_lim[level]; c--){
        Var x = var(trail[c]);
        variableDatabase.setVar(x, l_Undef);
        solver.branchingHeuristicManager.handleEventLitUnassigned(
            trail[c],
            solver.conflicts,
            c > trail_lim.last()
        );
    }

    // Decrease the size of the trail
    const int qhead = trail_lim[level];
    trail.shrink(trail.size() - trail_lim[level]);
    trail_lim.shrink(trail_lim.size() - level);

    // Add the assignments at 'level' to the queue
    solver.propagationQueue.batchEnqueue(trail, qhead);
}