#include <core/UnitPropagator.h>
#include <core/Solver.h>

namespace Minisat {

    UnitPropagator::UnitPropagator(Solver* s)
        : watches(s->ca)
        , propagation_budget(-1)
        , propagations(0)
        , variableDatabase(s->variableDatabase)
        , assignmentTrail(s->assignmentTrail)
        , propagationQueue(s->propagationQueue)
        , ca(s->ca)
        , solver(s)
    {}

    UnitPropagator::~UnitPropagator() {
        solver = nullptr;
    }

    void UnitPropagator::relocAll(ClauseAllocator& to) {
        // Clean watchers
        watches.cleanAll();

        // Reloc clauses in watchers
        for (int v = 0; v < variableDatabase.nVars(); v++) {
            Lit p = mkLit(v);
            relocWatchers(watches[ p], to);
            relocWatchers(watches[~p], to);
        }
    }

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
            if (variableDatabase.value(blocker) == l_True) {
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
            if (first != blocker && variableDatabase.value(first) == l_True) {
                *j++ = w;
                continue;
            }

            // Look for new watch:
            for (int k = 2; k < c.size(); k++) {
                if (variableDatabase.value(c[k]) != l_False) {
                    std::swap(c[1], c[k]);
                    watches[~c[1]].push(w);
                    goto NextClause;
                }
            }

            // Did not find watch -- clause is unit under assignment:
            *j++ = w;
            if (
                variableDatabase.value(first) == l_False ||
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
        solver->simpDB_props -= num_props;

        return confl;
    }

    void UnitPropagator::enforceWatcherInvariant(CRef cr, int i_undef, int i_max) {
        // Move unassigned literal to c[0]
        Clause& c = ca[cr];
        Lit x = c[i_undef], max = c[i_max];
        if (c.size() == 2) {
            // Don't need to touch watchers for binary clauses
            if (variableDatabase.value(c[0]) == l_False) std::swap(c[0], c[1]);
        } else {
            // Swap unassigned literal to index 0 and highest-level literal to index 1,
            // replacing watchers as necessary
            Lit c0 = c[0], c1 = c[1];
            if (i_max == 0 || i_undef == 1) std::swap(c[0], c[1]);

            if (i_max > 1) {
                remove(watches[~c[1]], Watcher(cr, c0));
                std::swap(c[1], c[i_max]);
                watches[~max].push(Watcher(cr, x));
            }

            if (i_undef > 1) {
                remove(watches[~c[0]], Watcher(cr, c1));
                std::swap(c[0], c[i_undef]);
                watches[~x].push(Watcher(cr, max));
            }
        }
    }
}