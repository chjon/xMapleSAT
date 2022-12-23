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
    solver->branchingHeuristicManager.handleEventLitAssigned(p, solver->conflicts);
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
            variableDatabase.setVar(x, l_Undef);
            solver->branchingHeuristicManager.handleEventLitUnassigned(trail[c], solver->conflicts, c > trail_lim.last());
        }

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