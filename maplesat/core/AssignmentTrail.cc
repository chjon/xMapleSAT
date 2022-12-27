#include "core/AssignmentTrail.h"
#include "core/Solver.h"

using namespace Minisat;

AssignmentTrail::AssignmentTrail(Solver* s)
    : variableDatabase(s->variableDatabase)
    , ca(s->ca)
    , solver(s)
{}

void AssignmentTrail::assign(Lit p, CRef from) {
    assert(variableDatabase.value(p) == l_Undef);
    const Var v = var(p);
    solver->branchingHeuristicManager.handleEventLitAssigned(p, solver->conflicts);
    variableDatabase.setVar(v, lbool(!sign(p)));
    vardata[v] = mkVarData(from, decisionLevel());
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
        for (int c = trail.size() - 1; c >= trail_lim[level]; c--){
            Var      x  = var(trail[c]);
            variableDatabase.setVar(x, l_Undef);
            solver->branchingHeuristicManager.handleEventLitUnassigned(trail[c], solver->conflicts, c > trail_lim.last());
        }

        const int qhead = trail_lim[level];
        trail.shrink(trail.size() - trail_lim[level]);
        trail_lim.shrink(trail_lim.size() - level);

        // Add the assignments at 'level' to the queue
        solver->propagationQueue.batchEnqueue(trail, qhead);
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