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

void UnitPropagator::relocAll(ClauseAllocator& to) {
    // Clean watchers
    watches_bin.cleanAll();
    watches    .cleanAll();

    // Reloc clauses in watchers
    for (int v = 0; v < assignmentTrail.nVars(); v++) {
        Lit p = mkLit(v);
        relocWatchers(watches_bin[ p], to);
        relocWatchers(watches_bin[~p], to);
        relocWatchers(watches    [ p], to);
        relocWatchers(watches    [~p], to);
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
        if (confl != CRef_Undef) {
        #ifndef LOOSE_PROP_STAT
            break;
        #else
            return confl;
        #endif
        }
    }

    propagations += num_props;
    solver.simpDB_props -= num_props;

    return confl;
}

CRef UnitPropagator::simplePropagate() {
    CRef    confl = CRef_Undef;
    int     num_props = 0;
    watches.cleanAll();
    watches_bin.cleanAll();

    Lit p; // 'p' is enqueued fact to propagate.
    while ((p = propagationQueue.getNext()) != lit_Undef) {
        vec<Watcher>&  ws = watches[p];
        Watcher        *i, *j, *end;
        num_props++;

        // First, Propagate binary clauses
        vec<Watcher>&  wbin = watches_bin[p];

        for (int k = 0; k < wbin.size(); k++) {
            Lit imp = wbin[k].blocker;
            if (assignmentTrail.value(imp) == l_False) {
                return wbin[k].cref;
            }

            if (assignmentTrail.value(imp) == l_Undef) {
                assignmentTrail.simpleUncheckEnqueue(imp, wbin[k].cref);
            }
        }

        for (i = j = (Watcher*)ws, end = i + ws.size(); i != end;) {
            // Try to avoid inspecting the clause:
            Lit blocker = i->blocker;
            if (assignmentTrail.value(blocker) == l_True) {
                *j++ = *i++; continue;
            }

            // Make sure the false literal is data[1]:
            CRef     cr = i->cref;
            Clause&  c = ca[cr];
            Lit      false_lit = ~p;
            if (c[0] == false_lit)
                c[0] = c[1], c[1] = false_lit;
            assert(c[1] == false_lit);
            //  i++;

            // If 0th watch is true, then clause is already satisfied.
            // However, 0th watch is not the blocker, make it blocker using a new watcher w
            // why not simply do i->blocker=first in this case?
            Lit first = c[0];
            //  Watcher w     = Watcher(cr, first);
            if (first != blocker && assignmentTrail.value(first) == l_True) {
                i->blocker = first;
                *j++ = *i++; continue;
            } else {  // ----------------- DEFAULT  MODE (NOT INCREMENTAL)
                for (int k = 2; k < c.size(); k++) {
                    if (assignmentTrail.value(c[k]) != l_False) {
                        // watcher i is abandonned using i++, because cr watches now ~c[k] instead of p
                        // the blocker is first in the watcher. However,
                        // the blocker in the corresponding watcher in ~first is not c[1]
                        Watcher w = Watcher(cr, first); i++;
                        c[1] = c[k]; c[k] = false_lit;
                        watches[~c[1]].push(w);
                        goto NextClause;
                    }
                }
            }

            // Did not find watch -- clause is unit under assignment:
            i->blocker = first;
            *j++ = *i++;
            if (assignmentTrail.value(first) == l_False) {
                confl = cr;
                propagationQueue.clear();
                break;
            } else {
                assignmentTrail.simpleUncheckEnqueue(first, cr);
            }
NextClause:;
        }
        // Copy the remaining watches:
        while (i < end) *j++ = *i++;
        ws.shrink(i - j);
    }

    return confl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

inline CRef UnitPropagator::propagateSingleBinary(Lit p) {
    vec<Watcher>& ws_bin = watches_bin[p];
    for (int k = 0; k < ws_bin.size(); k++){
        Lit the_other = ws_bin[k].blocker;
        if (assignmentTrail.value(the_other) == l_False){
            return ws_bin[k].cref;
        } else if (assignmentTrail.value(the_other) == l_Undef) {
            propagationQueue.enqueue(the_other, ws_bin[k].cref);
        }
    }

    return CRef_Undef;
}

inline CRef UnitPropagator::propagateSingleNonBinary(Lit p) {
    // Iterate through the watches for p, using pointers instead of array indices for speed
    CRef confl = CRef_Undef;
    Watcher *i, *j, *end;
    vec<Watcher>& ws = watches[p];

    for (i = j = (Watcher*)ws, end = i + ws.size(); i != end;) {
        // Try to avoid inspecting the clause:
        Lit blocker = i->blocker;
        if (assignmentTrail.value(blocker) == l_True) {
            *j++ = *i++; continue;
        }

        // Make sure the false literal is data[1]:
        CRef     cr        = i->cref;
        Clause&  c         = ca[cr];
        Lit      false_lit = ~p;
        if (c[0] == false_lit)
            c[0] = c[1], c[1] = false_lit;
        assert(c[1] == false_lit);
        i++;

        // If 0th watch is true, then clause is already satisfied.
        Lit     first = c[0];
        Watcher w = Watcher(cr, first);
        if (first != blocker && assignmentTrail.value(first) == l_True){
            *j++ = w; continue;
        }

        // Look for new watch:
        for (int k = 2; k < c.size(); k++) {
            if (assignmentTrail.value(c[k]) != l_False) {
                c[1] = c[k]; c[k] = false_lit;
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
    while (i < end) *j++ = *i++;
    ws.shrink(i - j);
    return confl;
}

inline CRef UnitPropagator::propagateSingle(Lit p) {
    // Propagate binary clauses first.
    CRef confl = propagateSingleBinary(p);
    if (confl != CRef_Undef) return confl;

    // Propagate non-binary clauses second.
    return propagateSingleNonBinary(p);
}