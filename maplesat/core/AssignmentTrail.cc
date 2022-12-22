#include "core/AssignmentTrail.h"
#include "core/Solver.h"

using namespace Minisat;

AssignmentTrail::AssignmentTrail(Solver* s)
    : variableDatabase(s->variableDatabase)
    , ca(s->ca)
    , solver(s)
{}

void AssignmentTrail::newVar(Var v) {
    vardata.push(mkVarData(CRef_Undef, 0));
    trail  .capacity(v + 1);
}

void AssignmentTrail::uncheckedEnqueue(Lit p, CRef from) {
    assert(variableDatabase.value(p) == l_Undef);
    solver->picked[var(p)] = solver->conflicts;
#if ANTI_EXPLORATION
    uint64_t age = solver->conflicts - solver->canceled[var(p)];
    if (age > 0) {
        double decay = pow(0.95, age);
        solver->activity[var(p)] *= decay;
#if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
        if (solver->order_heap_extlvl.inHeap(var(p)))
            solver->order_heap_extlvl.increase(var(p));
        if (solver->order_heap_degree.inHeap(var(p)))
            solver->order_heap_degree.increase(var(p));
#else
        if (solver->order_heap.inHeap(var(p))) {
            solver->order_heap.increase(var(p));
        }
#endif
    }
#endif
    solver->conflicted[var(p)] = 0;
#if ALMOST_CONFLICT
    solver->almost_conflicted[var(p)] = 0;
#endif
    variableDatabase.setVar(var(p), lbool(!sign(p)));
    vardata[var(p)] = mkVarData(from, decisionLevel());
    trail.push_(p);
}

void AssignmentTrail::relocAll(ClauseAllocator& to) {
    for (int i = 0; i < trail.size(); i++) {
        Var v = var(trail[i]);
        if (reason(v) != CRef_Undef && (ca[reason(v)].reloced() || locked(ca[reason(v)])))
            ca.reloc(vardata[v].reason, to);
    }
}

// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void AssignmentTrail::cancelUntil(int level) {
    if (decisionLevel() > level){
        for (int c = trail.size()-1; c >= trail_lim[level]; c--){
            Var      x  = var(trail[c]);
            uint64_t age = solver->conflicts - solver->picked[x];
            if (age > 0) {
                double reward = ((double) solver->conflicted[x]) / ((double) age);
#if BRANCHING_HEURISTIC == LRB
#if ALMOST_CONFLICT
                double adjusted_reward = ((double) (solver->conflicted[x] + solver->almost_conflicted[x])) / ((double) age);
#else
                double adjusted_reward = reward;
#endif
                double old_activity = solver->activity[x];
                solver->activity[x] = solver->step_size * adjusted_reward + ((1 - solver->step_size) * old_activity);
#if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
                auto& order_heap = solver->extensionLevel[x] ? solver->order_heap_extlvl : solver->order_heap_degree;
#else
                auto& order_heap = solver->order_heap;
#endif
                if (order_heap.inHeap(x)) {
                    if (solver->activity[x] > old_activity)
                        order_heap.decrease(x);
                    else
                        order_heap.increase(x);
                }
#endif
                solver->total_actual_rewards[x] += reward;
                solver->total_actual_count[x] ++;
            }
#if ANTI_EXPLORATION
            solver->canceled[x] = solver->conflicts;
#endif
            variableDatabase.setVar(x, l_Undef);
            if (solver->phase_saving > 1 || (solver->phase_saving == 1) && c > trail_lim.last())
                solver->polarity[x] = sign(trail[c]);
            solver->insertVarOrder(x); }

        solver->unitPropagator.setQueueHead(trail_lim[level]);
        trail.shrink(trail.size() - trail_lim[level]);
        trail_lim.shrink(trail_lim.size() - level);
    }
}

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