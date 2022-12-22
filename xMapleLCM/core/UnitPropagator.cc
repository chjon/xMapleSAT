#include <core/UnitPropagator.h>
#include <core/Solver.h>

namespace Minisat {

    UnitPropagator::UnitPropagator(Solver* s)
        : watches_bin(s->ca)
        , watches(s->ca)
        , bcp_order_heap(s->branchingHeuristicManager.getActivityVSIDS())
        , qhead(0)
        , propagation_budget(-1)
        , propagations(0)
        , s_propagations(0)
        , variableDatabase(s->variableDatabase)
        , ca(s->ca)
        , solver(s)
    {}

    UnitPropagator::~UnitPropagator() {
        solver = nullptr;
    }

    void UnitPropagator::relocAll(ClauseAllocator& to) {
        // Clean watchers
        watches_bin.cleanAll();
        watches.cleanAll();

        // Reloc clauses in watchers
        for (int v = 0; v < variableDatabase.nVars(); v++) {
            Lit p = mkLit(v);
            relocWatchers(watches_bin[ p], to);
            relocWatchers(watches_bin[~p], to);
            relocWatchers(watches    [ p], to);
            relocWatchers(watches    [~p], to);
        }
    }

    inline CRef UnitPropagator::propagate_single(Lit p) {
        CRef confl = CRef_Undef;
        vec<Watcher>&  ws  = watches[p];
        Watcher        *i, *j, *end;

        vec<Watcher>& ws_bin = watches_bin[p];  // Propagate binary clauses first.
        for (int k = 0; k < ws_bin.size(); k++){
            Lit the_other = ws_bin[k].blocker;
            if (variableDatabase.value(the_other) == l_True) {
                continue;
    #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
            } else if (variableDatabase.value(the_other) == l_False || bcpValue(the_other) == l_False){
                // Found a conflict!
                // Make sure conflicting literal is on the trail
                if (variableDatabase.value(the_other) == l_Undef)
                    solver->uncheckedEnqueue(~the_other, solver->vardata[var(the_other)].reason);

                // Clear the propagation queue
                for (int k = 0; k < bcp_order_heap.size(); k++)
                    bcp_assigns[bcp_order_heap[k] >> 1] = l_Undef;
                bcp_order_heap.clear();
    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
            } else if (variableDatabase.value(the_other) == l_False) {
                bcp_order_heap.clear();
    #else
            } else if (variableDatabase.value(the_other) == l_False) {
    #endif
                return ws_bin[k].cref;
            } else {
    #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
                // Queue literal for propagation
                if (bcpValue(the_other) == l_Undef) {
                    bcp_assigns[var(the_other)] = lbool(!sign(the_other));
                    solver->vardata[var(the_other)].reason = ws_bin[k].cref;
                    bcp_order_heap.insert(the_other.x);
                }
    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
                bcp_order_heap.insert(the_other.x);
                solver->uncheckedEnqueue(the_other, ws_bin[k].cref);
    #else
                solver->uncheckedEnqueue(the_other, ws_bin[k].cref);
    #endif
            }
        }

        for (i = j = (Watcher*)ws, end = i + ws.size();  i != end;){
            // Try to avoid inspecting the clause:
            Lit blocker = i->blocker;
            if (variableDatabase.value(blocker) == l_True){
                *j++ = *i++; continue; }

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
            Watcher w     = Watcher(cr, first);
            if (first != blocker && variableDatabase.value(first) == l_True){
                *j++ = w; continue; }

            // Look for new watch:
            for (int k = 2; k < c.size(); k++)
                if (variableDatabase.value(c[k]) != l_False){
                    c[1] = c[k]; c[k] = false_lit;
                    watches[~c[1]].push(w);
                    goto NextClause; }

            // Did not find watch -- clause is unit under assignment:
            *j++ = w;
    #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
            if (variableDatabase.value(first) == l_False || bcpValue(first) == l_False){
                // Found a conflict!
                // Make sure conflicting literal is on the trail
                if (variableDatabase.value(first) == l_Undef)
                    solver->uncheckedEnqueue(~first, solver->vardata[var(first)].reason);

                // Clear the propagation queue
                for (int k = 0; k < bcp_order_heap.size(); k++) {
                    Lit p; p.x = bcp_order_heap[k];
                    bcp_assigns[var(p)] = l_Undef;
                }
                bcp_order_heap.clear();
    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
            if (variableDatabase.value(first) == l_False) {
                bcp_order_heap.clear();
    #else
            if (variableDatabase.value(first) == l_False){
    #endif
                confl = cr;
                // Copy the remaining watches:
                while (i < end)
                    *j++ = *i++;
            } else {
    #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
                // Queue literal for propagation
                if (bcpValue(first) == l_Undef) {
                    bcp_assigns[var(first)] = lbool(!sign(first));
                    solver->vardata[var(first)].reason = cr;
                    bcp_order_heap.insert(first.x);
                }
    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
                bcp_order_heap.insert(first.x);
                solver->uncheckedEnqueue(first, cr);
    #else
                solver->uncheckedEnqueue(first, cr);
    #endif
            }

    NextClause:;
        }
        ws.shrink(i - j);
        return confl;
    }

    /*_________________________________________________________________________________________________
    |
    |  propagate : [void]  ->  [Clause*]
    |  
    |  Description:
    |    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
    |    otherwise CRef_Undef.
    |  
    |    Post-conditions:
    |      * the propagation queue is empty, even if there was a conflict.
    |________________________________________________________________________________________________@*/
    CRef UnitPropagator::propagate() {
        CRef    confl     = CRef_Undef;
        int     num_props = 0;
        Lit     p         = lit_Undef;
        watches.cleanAll();
        watches_bin.cleanAll();

    #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
        while (qhead < solver->trail.size() || !bcp_order_heap.empty()){
            if (qhead == solver->trail.size()) {
                p.x = bcp_order_heap.removeMin();
                solver->uncheckedEnqueue(p, solver->vardata[var(p)].reason);
                bcp_assigns[var(p)] = l_Undef;
            }
            p = solver->trail[qhead++];
    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
        for (int i = qhead; i < solver->trail.size(); i++)
            bcp_order_heap.insert(solver->trail[i].x);

        while (!bcp_order_heap.empty()) {
            p.x = bcp_order_heap.removeMin();
    #else
        while (qhead < solver->trail.size()){
            p = solver->trail[qhead++];
    #endif
            num_props++;
            confl = propagate_single(p);
            if (confl != CRef_Undef) {
                break;
            }
        }

        propagations += num_props;
        solver->simpDB_props -= num_props;
        qhead = solver->trail.size();

        return confl;
    }

    CRef UnitPropagator::simplePropagate() {
        CRef    confl = CRef_Undef;
        int     num_props = 0;
        watches.cleanAll();
        watches_bin.cleanAll();
        while (qhead < solver->trail.size())
        {
            Lit            p = solver->trail[qhead++];     // 'p' is enqueued fact to propagate.
            vec<Watcher>&  ws = watches[p];
            Watcher        *i, *j, *end;
            num_props++;


            // First, Propagate binary clauses
            vec<Watcher>&  wbin = watches_bin[p];

            for (int k = 0; k<wbin.size(); k++)
            {

                Lit imp = wbin[k].blocker;

                if (variableDatabase.value(imp) == l_False)
                {
                    return wbin[k].cref;
                }

                if (variableDatabase.value(imp) == l_Undef)
                {
                    solver->simpleUncheckEnqueue(imp, wbin[k].cref);
                }
            }
            for (i = j = (Watcher*)ws, end = i + ws.size(); i != end;)
            {
                // Try to avoid inspecting the clause:
                Lit blocker = i->blocker;
                if (variableDatabase.value(blocker) == l_True)
                {
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
                Lit     first = c[0];
                //  Watcher w     = Watcher(cr, first);
                if (first != blocker && variableDatabase.value(first) == l_True)
                {
                    i->blocker = first;
                    *j++ = *i++; continue;
                }
                else
                {  // ----------------- DEFAULT  MODE (NOT INCREMENTAL)
                    for (int k = 2; k < c.size(); k++)
                    {

                        if (variableDatabase.value(c[k]) != l_False)
                        {
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
                if (variableDatabase.value(first) == l_False)
                {
                    confl = cr;
                    qhead = solver->trail.size();
                    // Copy the remaining watches:
                    while (i < end)
                        *j++ = *i++;
                }
                else
                {
                    solver->simpleUncheckEnqueue(first, cr);
                }
    NextClause:;
            }
            ws.shrink(i - j);
        }

        s_propagations += num_props;

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