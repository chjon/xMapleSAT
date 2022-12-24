#include <core/UnitPropagator.h>
#include <core/Solver.h>

namespace Minisat {

    UnitPropagator::UnitPropagator(Solver* s)
        : watches(s->ca)
#if BCP_PRIORITY_MODE != BCP_PRIORITY_IMMEDIATE
        , bcp_order_heap(s->activity)
#endif
        , qhead(0)
        , propagation_budget(-1)
        , propagations(0)
        , variableDatabase(s->variableDatabase)
        , assignmentTrail(s->assignmentTrail)
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

    // inline CRef UnitPropagator::propagate_single(Lit p) {
        // CRef confl = CRef_Undef;
    //     vec<Watcher>&  ws  = watches[p];
    //     Watcher        *i, *j, *end;

    //     vec<Watcher>& ws_bin = watches_bin[p];  // Propagate binary clauses first.
    //     for (int k = 0; k < ws_bin.size(); k++){
    //         Lit the_other = ws_bin[k].blocker;
    //         if (variableDatabase.value(the_other) == l_True) {
    //             continue;
    // #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
    //         } else if (variableDatabase.value(the_other) == l_False || bcpValue(the_other) == l_False){
    //             // Found a conflict!
    //             // Make sure conflicting literal is on the trail
    //             if (variableDatabase.value(the_other) == l_Undef)
    //                 solver->uncheckedEnqueue(~the_other, solver->vardata[var(the_other)].reason);

    //             // Clear the propagation queue
    //             for (int k = 0; k < bcp_order_heap.size(); k++)
    //                 bcp_assigns[bcp_order_heap[k] >> 1] = l_Undef;
    //             bcp_order_heap.clear();
    // #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
    //         } else if (variableDatabase.value(the_other) == l_False) {
    //             bcp_order_heap.clear();
    // #else
    //         } else if (variableDatabase.value(the_other) == l_False) {
    // #endif
    //             return ws_bin[k].cref;
    //         } else {
    // #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
    //             // Queue literal for propagation
    //             if (bcpValue(the_other) == l_Undef) {
    //                 bcp_assigns[var(the_other)] = lbool(!sign(the_other));
    //                 solver->vardata[var(the_other)].reason = ws_bin[k].cref;
    //                 bcp_order_heap.insert(the_other.x);
    //             }
    // #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
    //             bcp_order_heap.insert(the_other.x);
    //             solver->uncheckedEnqueue(the_other, ws_bin[k].cref);
    // #else
    //             solver->uncheckedEnqueue(the_other, ws_bin[k].cref);
    // #endif
    //         }
    //     }

    //     for (i = j = (Watcher*)ws, end = i + ws.size();  i != end;){
    //         // Try to avoid inspecting the clause:
    //         Lit blocker = i->blocker;
    //         if (variableDatabase.value(blocker) == l_True){
    //             *j++ = *i++; continue; }

    //         // Make sure the false literal is data[1]:
    //         CRef     cr        = i->cref;
    //         Clause&  c         = ca[cr];
    //         Lit      false_lit = ~p;
    //         if (c[0] == false_lit)
    //             c[0] = c[1], c[1] = false_lit;
    //         assert(c[1] == false_lit);
    //         i++;

    //         // If 0th watch is true, then clause is already satisfied.
    //         Lit     first = c[0];
    //         Watcher w     = Watcher(cr, first);
    //         if (first != blocker && variableDatabase.value(first) == l_True){
    //             *j++ = w; continue; }

    //         // Look for new watch:
    //         for (int k = 2; k < c.size(); k++)
    //             if (variableDatabase.value(c[k]) != l_False){
    //                 c[1] = c[k]; c[k] = false_lit;
    //                 watches[~c[1]].push(w);
    //                 goto NextClause; }

    //         // Did not find watch -- clause is unit under assignment:
    //         *j++ = w;
    // #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
    //         if (variableDatabase.value(first) == l_False || bcpValue(first) == l_False){
    //             // Found a conflict!
    //             // Make sure conflicting literal is on the trail
    //             if (variableDatabase.value(first) == l_Undef)
    //                 solver->uncheckedEnqueue(~first, solver->vardata[var(first)].reason);

    //             // Clear the propagation queue
    //             for (int k = 0; k < bcp_order_heap.size(); k++) {
    //                 Lit p; p.x = bcp_order_heap[k];
    //                 bcp_assigns[var(p)] = l_Undef;
    //             }
    //             bcp_order_heap.clear();
    // #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
    //         if (variableDatabase.value(first) == l_False) {
    //             bcp_order_heap.clear();
    // #else
    //         if (variableDatabase.value(first) == l_False){
    // #endif
    //             confl = cr;
    //             // Copy the remaining watches:
    //             while (i < end)
    //                 *j++ = *i++;
    //         } else {
    // #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
    //             // Queue literal for propagation
    //             if (bcpValue(first) == l_Undef) {
    //                 bcp_assigns[var(first)] = lbool(!sign(first));
    //                 solver->vardata[var(first)].reason = cr;
    //                 bcp_order_heap.insert(first.x);
    //             }
    // #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
    //             bcp_order_heap.insert(first.x);
    //             solver->uncheckedEnqueue(first, cr);
    // #else
    //             solver->uncheckedEnqueue(first, cr);
    // #endif
    //         }

    // NextClause:;
    //     }
    //     ws.shrink(i - j);
    //     return confl;
    // }

    // /*_________________________________________________________________________________________________
    // |
    // |  propagate : [void]  ->  [Clause*]
    // |  
    // |  Description:
    // |    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
    // |    otherwise CRef_Undef.
    // |  
    // |    Post-conditions:
    // |      * the propagation queue is empty, even if there was a conflict.
    // |________________________________________________________________________________________________@*/
    // CRef UnitPropagator::propagate() {
    //     CRef    confl     = CRef_Undef;
    //     int     num_props = 0;
    //     Lit     p         = lit_Undef;
    //     watches.cleanAll();

    // #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
    //     while (qhead < solver->trail.size() || !bcp_order_heap.empty()){
    //         if (qhead == solver->trail.size()) {
    //             p.x = bcp_order_heap.removeMin();
    //             solver->uncheckedEnqueue(p, solver->vardata[var(p)].reason);
    //             bcp_assigns[var(p)] = l_Undef;
    //         }
    //         p = solver->trail[qhead++];
    // #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
    //     for (int i = qhead; i < solver->trail.size(); i++)
    //         bcp_order_heap.insert(solver->trail[i].x);

    //     while (!bcp_order_heap.empty()) {
    //         p.x = bcp_order_heap.removeMin();
    // #else
    //     while (qhead < solver->trail.size()){
    //         p = solver->trail[qhead++];
    // #endif
    //         num_props++;
    //         confl = propagate_single(p);
    //         if (confl != CRef_Undef) {
    //             break;
    //         }
    //     }

    //     propagations += num_props;
    //     solver->simpDB_props -= num_props;
    //     qhead = solver->trail.size();

        // return confl;
    // }

    inline void UnitPropagator::enqueue(Lit p, CRef from) {
    #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
        // Queue literal for propagation
        if (bcpValue(first) == l_Undef) {
            bcp_assigns[var(first)] = lbool(!sign(first));
            solver->vardata[var(first)].reason = cr;
            bcp_order_heap.insert(first.x);
        }
    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
        bcp_order_heap.insert(first.x);
        assignmentTrail.uncheckedEnqueue(p, from);
    #else
        assignmentTrail.uncheckedEnqueue(p, from);
    #endif
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
            if (variableDatabase.value(first) == l_False) {
                // All literals falsified!
                confl = cr;
                break;
            } else {
                // Propagate literal
                enqueue(first, cr);
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
        
        while (qhead < assignmentTrail.nAssigns()) {
            num_props++;
            Lit p = assignmentTrail[qhead++]; // 'p' is enqueued fact to propagate.
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