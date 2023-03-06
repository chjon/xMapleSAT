#include "core/Solver.h"

using namespace Minisat;
using VarData = Solver::VarData;
using VarOrderLt = Solver::VarOrderLt;

static inline void updateCHB(
    Var x, uint64_t conflicts
    , vec<double>& activity_CHB
    , Heap<VarOrderLt>& order_heap_CHB
    , vec<uint32_t>& picked
    , vec<uint32_t>& conflicted
    , vec<uint32_t>& almost_conflicted
#ifdef ANTI_EXPLORATION
    , vec<uint32_t>& canceled
#endif
) {
    picked[x] = conflicts;
    conflicted[x] = 0;
    almost_conflicted[x] = 0;
#ifdef ANTI_EXPLORATION
    const uint32_t age = conflicts - canceled[x];
    if (age > 0){
        double decay = pow(0.95, age);
        activity_CHB[x] *= decay;
        if (order_heap_CHB.inHeap(x))
            order_heap_CHB.increase(x);
    }
#endif
}

static inline void enqueueGreedy(
    Lit p, int level, CRef from,
    vec<lbool>& assigns, vec<VarData>& vardata, vec<Lit>& trail
) {
    const Var x = var(p);
    assigns[x] = lbool(!sign(p));
    vardata[x] = Solver::mkVarData(from, level);
    trail.push_(p);
}

static inline void enqueueOutOfOrder(
    Lit p, int level, CRef from,
    vec<lbool>& assigns, vec<VarData>& vardata, Heap<VarOrderLt>& bcp_heap, int qhead, vec<Lit>& trail
) {
    // Store reason clause and decision level
    const Var x = var(p);
    assigns[x] = lbool(!sign(p));
    vardata[x] = Solver::mkVarData(from, level);

    if (qhead == trail.size()) {
        bcp_heap.insert(x);
        qhead++;
    }
    trail.push_(p);
}

static inline void enqueueDelayed(
    Lit p, int level, CRef from,
    vec<lbool>& bcp_assigns, vec<VarData>& vardata, Heap<VarOrderLt>& bcp_heap
) {
    const Var x = var(p);
    if (bcp_assigns[x] == l_Undef) { // Needed in case bcp_assigns[x] == l_True
        bcp_assigns[x] = lbool(!sign(p));
        vardata[x] = Solver::mkVarData(from, level);
        bcp_heap.insert(x);
    }
}

void Solver::uncheckedEnqueue(Lit p, CRef from) {
    assert(value(p) == l_Undef);

    // Update branching heuristics
    if (!VSIDS) updateCHB(
        var(p), conflicts, activity_CHB, order_heap_CHB, picked, conflicted, almost_conflicted
    #ifdef ANTI_EXPLORATION
        , canceled
    #endif
    );

    // Add literal to propagation queue
    const int level = decisionLevel();
    switch (bcp_mode) {
        case BCPMode::DELAYED:
            enqueueDelayed(p, level, from, bcp_assigns, vardata, bcp_heap); break;
        case BCPMode::OUT_OF_ORDER:
            enqueueOutOfOrder(p, level, from, assigns, vardata, bcp_heap, qhead, trail); break;
        default:
            enqueueGreedy(p, level, from, assigns, vardata, trail); break;
    }
}

static inline Lit getNextGreedy(
    int& qhead,
    vec<Lit>& trail
) {
    return qhead < trail.size() ? trail[qhead++] : lit_Undef;
}

static inline Lit getNextOutOfOrder(
    vec<lbool>& assigns,
    Heap<VarOrderLt>& bcp_heap,
    int& qhead,
    vec<Lit>& trail
) {
    // Check if the heap is empty
    if (bcp_heap.size() == 0) return getNextGreedy(qhead, trail);

    // Propagate from the heap
    const Var x = bcp_heap.removeMin();
    return mkLit(x, assigns[x] == l_False);
}

static inline Lit getNextDelayed(
    vec<lbool>& assigns,
    vec<lbool>& bcp_assigns,
    Heap<VarOrderLt>& bcp_heap,
    int& qhead,
    vec<Lit>& trail
) {
    // Prioritize propagating from the trail first
    if (qhead < trail.size()) return trail[qhead++];

    // Check if the heap is empty
    if (bcp_heap.size() == 0) return lit_Undef;
    
    // Propagate from the heap
    qhead++;
    const Var x = bcp_heap.removeMin();
    assigns[x] = bcp_assigns[x];
    bcp_assigns[x] = l_Undef;
    const Lit p = mkLit(x, assigns[x] == l_False);
    trail.push_(p);
    return p;
}

static inline Lit getNext(
    BCPMode bcp_mode,
    vec<lbool>& assigns,
    vec<lbool>& bcp_assigns,
    Heap<VarOrderLt>& bcp_heap,
    int& qhead,
    vec<Lit>& trail
) {
    switch (bcp_mode) {
        case BCPMode::DELAYED:
            return getNextDelayed(assigns, bcp_assigns, bcp_heap, qhead, trail);
        case BCPMode::OUT_OF_ORDER:
            return getNextOutOfOrder(assigns, bcp_heap, qhead, trail);
        default:
            return getNextGreedy(qhead, trail);
    }
}

static inline bool causesConflictGreedy(Lit p, vec<lbool>& assigns) {
    return (assigns[var(p)] ^ sign(p)) == l_False;
}

static inline bool causesConflictDelayed(
    Lit p,
    vec<lbool>& assigns,
    vec<lbool>& bcp_assigns,
    vec<Lit>& trail
) {
    const Var x = var(p);
    if ((assigns[x] ^ sign(p)) == l_False) return true;
    if ((bcp_assigns[x] ^ sign(p)) != l_False) return false;

    assigns[x] = bcp_assigns[x];
    bcp_assigns[x] = l_Undef;
    trail.push_(~p);
    return true;
}

static inline bool causesConflict(
    BCPMode bcp_mode,
    Lit p,
    vec<lbool>& assigns,
    vec<lbool>& bcp_assigns,
    vec<Lit>& trail
) {
    switch (bcp_mode) {
        case BCPMode::DELAYED:
            return causesConflictDelayed(p, assigns, bcp_assigns, trail);
        default:
            return causesConflictGreedy(p, assigns);
    }
}

static inline void clearQueueGreedy(
    int& qhead,
    vec<Lit>& trail
) {
    qhead = trail.size();
}

static inline void clearQueueOutOfOrder(
    Heap<VarOrderLt>& bcp_heap,
    int& qhead,
    vec<Lit>& trail
) {
    bcp_heap.clear();
    qhead = trail.size();
}

static inline void clearQueueDelayed(
    vec<lbool>& bcp_assigns,
    Heap<VarOrderLt>& bcp_heap,
    int& qhead,
    vec<Lit>& trail
) {
    for (int i = 0; i < bcp_heap.size(); i++)
        bcp_assigns[bcp_heap[i]] = l_Undef;
    bcp_heap.clear();
    qhead = trail.size();
}

static inline void clearQueue(
    BCPMode bcp_mode,
    vec<lbool>& bcp_assigns,
    Heap<VarOrderLt>& bcp_heap,
    int& qhead,
    vec<Lit>& trail
) {
    switch (bcp_mode) {
        case BCPMode::DELAYED:
            clearQueueDelayed(bcp_assigns, bcp_heap, qhead, trail); break;
        case BCPMode::OUT_OF_ORDER:
            clearQueueOutOfOrder(bcp_heap, qhead, trail); break;
        default:
            clearQueueGreedy(qhead, trail); break;
    }
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
CRef Solver::propagate()
{
    CRef    confl     = CRef_Undef;
    int     num_props = 0;
    watches.cleanAll();
    watches_bin.cleanAll();

    for (
        Lit p = getNext(bcp_mode, assigns, bcp_assigns, bcp_heap, qhead, trail);
        p != lit_Undef;
        p = getNext(bcp_mode, assigns, bcp_assigns, bcp_heap, qhead, trail)
    ) {
        vec<Watcher>&  ws  = watches[p];
        Watcher        *i, *j, *end;
        num_props++;

        vec<Watcher>& ws_bin = watches_bin[p];  // Propagate binary clauses first.
        for (int k = 0; k < ws_bin.size(); k++){
            Lit the_other = ws_bin[k].blocker;
            if (causesConflict(bcp_mode, the_other, assigns, bcp_assigns, trail)){
                clearQueue(bcp_mode, bcp_assigns, bcp_heap, qhead, trail);
                confl = ws_bin[k].cref;
#ifdef LOOSE_PROP_STAT
                return confl;
#else
                goto ExitProp;
#endif
            }else if(value(the_other) == l_Undef)
                uncheckedEnqueue(the_other, ws_bin[k].cref);
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
            if (causesConflict(bcp_mode, first, assigns, bcp_assigns, trail)){
                clearQueue(bcp_mode, bcp_assigns, bcp_heap, qhead, trail);
                confl = cr;
                // Copy the remaining watches:
                while (i < end)
                    *j++ = *i++;
            }else
                uncheckedEnqueue(first, cr);

NextClause:;
        }
        ws.shrink(i - j);
    }

#ifndef LOOSE_PROP_STAT
ExitProp:;
#endif
    propagations += num_props;
    simpDB_props -= num_props;

    return confl;
}
