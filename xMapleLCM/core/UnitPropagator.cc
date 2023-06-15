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
    , watches_bin(s.ca)
    , watches(s.ca)
    , propagation_budget(-1)
    , propagations(0)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// UTILITY FUNCTIONS

void UnitPropagator::enforceWatcherInvariant(CRef cr, int i_undef, int i_max) {
    // Move unassigned literal to c[0]
    Clause& c = ca[cr];
    Lit x = c[i_undef], max = c[i_max];
    if (c.size() == 2) {
        // Don't need to touch watchers for binary clauses
        if (assignmentTrail.value(c[0]) == l_False) std::swap(c[0], c[1]);
    } else {
        // Swap unassigned literal to index 0 and highest-level literal to index 1,
        // replacing watchers as necessary
        OccLists<Lit, vec<Watcher>, WatcherDeleted>& ws = watches;

        Lit c0 = c[0], c1 = c[1];
        if (i_max == 0 || i_undef == 1) std::swap(c[0], c[1]);

        if (i_max > 1) {
            remove(ws[~c[1]], Watcher(cr, c0));
            std::swap(c[1], c[i_max]);
            ws[~max].push(Watcher(cr, x));
        }

        if (i_undef > 1) {
            remove(ws[~c[0]], Watcher(cr, c1));
            std::swap(c[0], c[i_undef]);
            ws[~x].push(Watcher(cr, max));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// STATE MODIFICATION

static inline void relocWatchers(vec<Watcher>& ws, ClauseAllocator& from, ClauseAllocator& to) {
    for (int i = 0; i < ws.size(); i++) from.reloc(ws[i].cref, to);
}

void UnitPropagator::relocAll(ClauseAllocator& to) {
    // Clean watchers
    watches_bin.cleanAll();
    watches    .cleanAll();

    // Reloc clauses in watchers
    for (int v = 0; v < assignmentTrail.nVars(); v++) {
        Lit p = mkLit(v);
        relocWatchers(watches_bin[ p], ca, to);
        relocWatchers(watches_bin[~p], ca, to);
        relocWatchers(watches    [ p], ca, to);
        relocWatchers(watches    [~p], ca, to);
    }
}

CRef UnitPropagator::propagate() {
    switch (propagationQueue.current_bcpmode) {
        default:
        case BCPMode::IMMEDIATE : return genericPropagate<BCPMode::IMMEDIATE , false>();
        case BCPMode::DELAYED   : return genericPropagate<BCPMode::DELAYED   , false>();
        case BCPMode::OUTOFORDER: return genericPropagate<BCPMode::OUTOFORDER, false>();
    }
}

CRef UnitPropagator::simplePropagate() {
    switch (propagationQueue.current_bcpmode) {
        default:
        case BCPMode::IMMEDIATE : return genericPropagate<BCPMode::IMMEDIATE , true>();
        case BCPMode::DELAYED   : return genericPropagate<BCPMode::DELAYED   , true>();
        case BCPMode::OUTOFORDER: return genericPropagate<BCPMode::OUTOFORDER, true>();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS FOR propagate()

template <BCPMode bcpmode, bool simple>
static inline bool enqueue(PropagationQueue& propagationQueue, Lit p, CRef from) {
    if (simple) return propagationQueue.simpleEnqueue<bcpmode>(p, from);
    else        return propagationQueue.enqueue      <bcpmode>(p, from);
}

template <BCPMode bcpmode, bool simple>
static inline CRef propagateSingleBinary(
    PropagationQueue& pq,
    vec<Watcher>& ws_bin,
    const AssignmentTrail& at
) {
    for (int k = 0; k < ws_bin.size(); k++) {
        Lit the_other = ws_bin[k].blocker;
        if (at.value(the_other) == l_True) continue;
        if (
            at.value(the_other) == l_False ||
            !enqueue<bcpmode, simple>(pq, the_other, ws_bin[k].cref)
        ) {
            // All literals falsified!
            pq.clear();
            return ws_bin[k].cref;
        }
    }

    return CRef_Undef;
}

static inline bool findNewWatch(
    OccLists<Lit, vec<Watcher>, WatcherDeleted>& watches,
    const AssignmentTrail& at,
    Clause& c,
    CRef cr
) {
    for (int k = 2; k < c.size(); k++) {
        if (at.value(c[k]) != l_False) {
            std::swap(c[1], c[k]);
            watches[~c[1]].push(Watcher(cr, c[0]));
            return true;
        }
    }
    return false;
}

template <BCPMode bcpmode, bool simple>
inline CRef UnitPropagator::propagateSingleNonBinary(Lit p) {
    // Iterate through the watches for p, using pointers instead of array indices for speed
    CRef confl = CRef_Undef;
    Watcher *i, *j, *end;
    vec<Watcher>& ws = watches[p];

    for (i = j = (Watcher*)ws, end = i + ws.size(); i != end;) {
        // Try to avoid inspecting the clause:
        Lit blocker = i->blocker;
        if (assignmentTrail.value(blocker) == l_True) {
            *j++ = *i++;
            continue;
        }

        // Make sure the false literal is data[1]:
        CRef cr = i->cref;
        Clause& c = ca[cr];
        if (c[0] == ~p) std::swap(c[0], c[1]);
        assert(c[1] == ~p);

        // If 0th watch is true, then clause is already satisfied.
        // If 0th watch is not already the blocker, make it the blocker
        Lit first = c[0];
        if (first != blocker && assignmentTrail.value(first) == l_True) {
            i->blocker = first;
            *j++ = *i++;
            continue;
        }
        
        // Look for new watch
        if (findNewWatch(watches, assignmentTrail, c, cr)) {
            i++;
            continue;
        }

        // Did not find watch -- clause is unit under assignment:
        i->blocker = first;
        *j++ = *i++;

        if (
            assignmentTrail.value(first) == l_False ||
            !enqueue<bcpmode, simple>(propagationQueue, first, cr)
        ) {
            // All literals falsified!
            confl = cr;
            propagationQueue.clear();
            break;
        }
    }

    // Copy the remaining watches:
    while (i < end) *j++ = *i++;
    ws.shrink(i - j);
    return confl;
}

template <BCPMode bcpmode, bool simple>
inline CRef UnitPropagator::genericPropagate() {
    // Batch update propagations stat to avoid cost of data non-locality
    CRef confl = CRef_Undef;
    int num_props = 0;
    watches.cleanAll();
    watches_bin.cleanAll();

    Lit p; // 'p' is enqueued fact to propagate.
    while ((p = propagationQueue.getNext<bcpmode, simple>()) != lit_Undef) {
        num_props++;

        // Propagate binary clauses first.
        confl = propagateSingleBinary<bcpmode, simple>(propagationQueue, watches_bin[p], assignmentTrail);
        if (confl != CRef_Undef) {
            if (simple) {
                break;
            } else {
            #ifdef LOOSE_PROP_STAT
                return confl;
            #else
                break;
            #endif
            }
        }

        // Propagate non-binary clauses second.
        confl = propagateSingleNonBinary<bcpmode, simple>(p);
        if (confl != CRef_Undef) break;
    }

    if (!simple) {
        propagations += num_props;
        solver.simpDB_props -= num_props;
    }

    return confl;
}
