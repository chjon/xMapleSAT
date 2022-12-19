#include <core/PropagationComponent.h>
#include <core/Solver.h>

namespace Minisat {

    PropagationComponent::PropagationComponent(Solver* s)
        : watches_bin(s->ca)
        , watches(s->ca)
        , bcp_order_heap(s->activity_VSIDS)
        , propagation_budget(-1)
        , propagations(0)
        , s_propagations(0)
        , ca(s->ca)
        , solver(s)
    {}

    PropagationComponent::~PropagationComponent() {
        solver = nullptr;
    }

#ifdef TESTING
    inline void  PropagationComponent::set_value(Var x, lbool v, int l) {
        auto it = test_value.find(x);
        if (it == test_value.end()) test_value.insert(std::make_pair(x, std::make_pair(v, l)));
        else it->second = std::make_pair(v, l);
    }

    inline lbool PropagationComponent::value(Var x) const { auto it = test_value.find(x); return (it == test_value.end()) ? (l_Undef) : (it->second.first ); }
    inline lbool PropagationComponent::value(Lit p) const { return value(var(p)) ^ sign(p); }
#else
    inline lbool PropagationComponent::value(Var x) const { return solver->value(x); }
    inline lbool PropagationComponent::value(Lit p) const { return solver->value(p); }
#endif

    void PropagationComponent::relocAll(ClauseAllocator& to) {
        // Clean watchers
        watches_bin.cleanAll();
        watches.cleanAll();

        // Reloc clauses in watchers
        for (int v = 0; v < solver->nVars(); v++) {
            Lit p = mkLit(v);
            relocWatchers(watches_bin[ p], to);
            relocWatchers(watches_bin[~p], to);
            relocWatchers(watches    [ p], to);
            relocWatchers(watches    [~p], to);
        }
    }

    inline CRef PropagationComponent::propagate_single(Lit p) {
        CRef confl = CRef_Undef;
        vec<Watcher>&  ws  = watches[p];
        Watcher        *i, *j, *end;

        vec<Watcher>& ws_bin = watches_bin[p];  // Propagate binary clauses first.
        for (int k = 0; k < ws_bin.size(); k++){
            Lit the_other = ws_bin[k].blocker;
            if (value(the_other) == l_True) {
                continue;
    #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
            } else if (value(the_other) == l_False || bcpValue(the_other) == l_False){
                // Found a conflict!
                // Make sure conflicting literal is on the trail
                if (value(the_other) == l_Undef)
                    solver->uncheckedEnqueue(~the_other, solver->vardata[var(the_other)].reason);

                // Clear the propagation queue
                for (int k = 0; k < bcp_order_heap.size(); k++)
                    bcp_assigns[bcp_order_heap[k] >> 1] = l_Undef;
                bcp_order_heap.clear();
    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
            } else if (value(the_other) == l_False) {
                bcp_order_heap.clear();
    #else
            } else if (value(the_other) == l_False) {
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
            if (value(blocker) == l_True){
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
            if (first != blocker && value(first) == l_True){
                *j++ = w; continue; }

            // Look for new watch:
            for (int k = 2; k < c.size(); k++)
                if (value(c[k]) != l_False){
                    c[1] = c[k]; c[k] = false_lit;
                    watches[~c[1]].push(w);
                    goto NextClause; }

            // Did not find watch -- clause is unit under assignment:
            *j++ = w;
    #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
            if (value(first) == l_False || bcpValue(first) == l_False){
                // Found a conflict!
                // Make sure conflicting literal is on the trail
                if (value(first) == l_Undef)
                    solver->uncheckedEnqueue(~first, solver->vardata[var(first)].reason);

                // Clear the propagation queue
                for (int k = 0; k < bcp_order_heap.size(); k++) {
                    Lit p; p.x = bcp_order_heap[k];
                    bcp_assigns[var(p)] = l_Undef;
                }
                bcp_order_heap.clear();
    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
            if (value(first) == l_False) {
                bcp_order_heap.clear();
    #else
            if (value(first) == l_False){
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
    CRef PropagationComponent::propagate() {
        CRef    confl     = CRef_Undef;
        int     num_props = 0;
        Lit     p         = lit_Undef;
        watches.cleanAll();
        watches_bin.cleanAll();

    #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
        while (solver->qhead < solver->trail.size() || !bcp_order_heap.empty()){
            if (solver->qhead == solver->trail.size()) {
                p.x = bcp_order_heap.removeMin();
                solver->uncheckedEnqueue(p, solver->vardata[var(p)].reason);
                bcp_assigns[var(p)] = l_Undef;
            }
            p = solver->trail[solver->qhead++];
    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
        for (int i = solver->qhead; i < solver->trail.size(); i++)
            bcp_order_heap.insert(solver->trail[i].x);

        while (!bcp_order_heap.empty()) {
            p.x = bcp_order_heap.removeMin();
    #else
        while (solver->qhead < solver->trail.size()){
            p = solver->trail[solver->qhead++];
    #endif
            num_props++;
            confl = propagate_single(p);
            if (confl != CRef_Undef) {
                break;
            }
        }

        propagations += num_props;
        solver->simpDB_props -= num_props;
        solver->qhead = solver->trail.size();

        return confl;
    }

    CRef PropagationComponent::simplePropagate() {
        CRef    confl = CRef_Undef;
        int     num_props = 0;
        watches.cleanAll();
        watches_bin.cleanAll();
        while (solver->qhead < solver->trail.size())
        {
            Lit            p = solver->trail[solver->qhead++];     // 'p' is enqueued fact to propagate.
            vec<Watcher>&  ws = watches[p];
            Watcher        *i, *j, *end;
            num_props++;


            // First, Propagate binary clauses
            vec<Watcher>&  wbin = watches_bin[p];

            for (int k = 0; k<wbin.size(); k++)
            {

                Lit imp = wbin[k].blocker;

                if (value(imp) == l_False)
                {
                    return wbin[k].cref;
                }

                if (value(imp) == l_Undef)
                {
                    solver->simpleUncheckEnqueue(imp, wbin[k].cref);
                }
            }
            for (i = j = (Watcher*)ws, end = i + ws.size(); i != end;)
            {
                // Try to avoid inspecting the clause:
                Lit blocker = i->blocker;
                if (value(blocker) == l_True)
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
                if (first != blocker && value(first) == l_True)
                {
                    i->blocker = first;
                    *j++ = *i++; continue;
                }
                else
                {  // ----------------- DEFAULT  MODE (NOT INCREMENTAL)
                    for (int k = 2; k < c.size(); k++)
                    {

                        if (value(c[k]) != l_False)
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
                if (value(first) == l_False)
                {
                    confl = cr;
                    solver->qhead = solver->trail.size();
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

    void PropagationComponent::enforceWatcherInvariant(CRef cr, int i_undef, int i_max) {
        // Move unassigned literal to c[0]
        Clause& c = ca[cr];
        Lit x = c[i_undef], max = c[i_max];
        if (c.size() == 2) {
            // Don't need to touch watchers for binary clauses
            if (value(c[0]) == l_False) std::swap(c[0], c[1]);
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