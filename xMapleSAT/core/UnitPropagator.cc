/*******************************************************************************[UnitPropagator.cc]
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

#include <core/UnitPropagator.h>
#include <core/Solver.h>

using namespace Minisat;

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

UnitPropagator::UnitPropagator(Solver& s)
    //////////////////////
    // Solver references
    : propagationQueue(s.propagationQueue)
    , assignmentTrail(s.assignmentTrail)
    , ca(s.ca)
    , solver(s)
    
    /////////////////////
    // Member variables
    , watches(s.ca)
    , propagation_budget(-1)
    , propagations(0)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// STATE MODIFICATION

void UnitPropagator::relocAll(ClauseAllocator& to) {
    // Clean watchers
    watches.cleanAll();

    // Reloc clauses in watchers
    for (int v = 0; v < assignmentTrail.nVars(); v++) {
        Lit p = mkLit(v);
        relocWatchers(watches[ p], to);
        relocWatchers(watches[~p], to);
    }
}

CRef UnitPropagator::propagate() {
    // Batch update propagations stat to avoid cost of data non-locality
    int num_props = 0;

    CRef confl = CRef_Undef;
    watches.cleanAll();
    
    Lit p; // 'p' is enqueued fact to propagate.
    while ((p = propagationQueue.getNext()) != lit_Undef) {
        num_props++;
        confl = propagateSingle(p);
        if (confl != CRef_Undef) break;
    }

    propagations += num_props;
    solver.simpDB_props -= num_props;

    return confl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

inline CRef UnitPropagator::propagateSingle(Lit p) {
    CRef confl = CRef_Undef;

    // Iterate through the watches for p, using pointers instead of array indices for speed
    vec<Watcher>& ws = watches[p];
    Watcher *i, *j, *end;
    i = j = (Watcher*)ws;
    end = i + ws.size();
    while (i != end) {
        // Try to avoid inspecting the clause:
        Lit blocker = i->blocker;
        if (assignmentTrail.value(blocker) == l_True) {
            *j++ = *i++;
            continue;
        }

        // Make sure the false literal is data[1]:
        CRef    cr = (i++)->cref;
        Clause& c = ca[cr];
        if (c[0] == ~p) std::swap(c[0], c[1]);
        assert(c[1] == ~p);

        // If 0th watch is true, then clause is already satisfied.
        const Lit first = c[0];
        Watcher w = Watcher(cr, first);
        if (first != blocker && assignmentTrail.value(first) == l_True) {
            *j++ = w;
            continue;
        }

        // Look for new watch:
        for (int k = 2; k < c.size(); k++) {
            if (assignmentTrail.value(c[k]) != l_False) {
                std::swap(c[1], c[k]);
                watches[~c[1]].push(w);
                goto NextClause;
            }
        }

        // Did not find watch -- clause is unit under assignment:
        *j++ = w;
        if (
            assignmentTrail.value(first) == l_False ||
            !propagationQueue.enqueue(first, cr)
        ) {
            // All literals falsified!
            confl = cr;
            propagationQueue.clear();
            break;
        }

    NextClause:;
    }

    // Copy the remaining watches:
    while (i != end) *j++ = *i++;
    ws.shrink(i - j);
    return confl;
}