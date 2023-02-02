/********************************************************************[BranchingHeuristicManager.cc]
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

#include "core/BranchingHeuristicManager.h"
#include "core/Solver.h"

using namespace Minisat;

///////////////////////////////////////////////////////////////////////////////////////////////////
// OPTIONS

static const char* _cat = "CORE";

static DoubleOption opt_var_decay     (_cat, "var-decay",   "The variable activity decay factor",            0.80,     DoubleRange(0, false, 1, false));
static DoubleOption opt_step_size     (_cat, "step-size",   "Initial step size",                             0.40,     DoubleRange(0, false, 1, false));
static DoubleOption opt_step_size_dec (_cat, "step-size-dec","Step size decrement",                          0.000001, DoubleRange(0, false, 1, false));
static DoubleOption opt_min_step_size (_cat, "min-step-size","Minimal step size",                            0.06,     DoubleRange(0, false, 1, false));
static BoolOption   opt_rnd_init_act  (_cat, "rnd-init",    "Randomize the initial activity", false);
static IntOption    opt_phase_saving  (_cat, "phase-saving", "Controls the level of phase saving (0=none, 1=limited, 2=full)", 2, IntRange(0, 2));
static IntOption    opt_VSIDS_props_limit (_cat, "VSIDS-lim",  "specifies the number of propagations after which the solver switches between LRB and VSIDS(in millions).", 30, IntRange(1, INT32_MAX));

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

BranchingHeuristicManager::BranchingHeuristicManager(Solver& s)
    //////////////////////////
    // Heuristic configuration

    : order_heap_CHB     (VarOrderLt<double>{activity_CHB, false})
    , order_heap_VSIDS   (VarOrderLt<double>{activity_VSIDS, false})
    , order_heap_distance(VarOrderLt<double>{activity_distance, false})

    // VSIDS
    , var_inc(1)
    , var_decay(opt_var_decay)
    , timer(5000)

    // CHB
    , step_size    (opt_step_size)
    , step_size_dec(opt_step_size_dec)
    , min_step_size(opt_min_step_size)

    // LRB
    , VSIDS(
        BRANCHING_HEURISTIC == BRANCHING_HEURISTIC_VSIDS
    )

    // Distance
    , DISTANCE(
        BRANCHING_HEURISTIC == BRANCHING_HEURISTIC_DYNAMIC ||
        BRANCHING_HEURISTIC == BRANCHING_HEURISTIC_DISTANCE
    )
    , var_iLevel_inc(1)
    , my_var_decay(0.6)

    // Mode switching
    , switchMode(false)
    , prevPropagations(0)

    // Random
    , rnd_init_act(opt_rnd_init_act)

    // Phase saving
    , phase_saving(static_cast<PhaseSavingLevel>(static_cast<int>(opt_phase_saving)))

    // Mode switching
    , VSIDS_props_limit(opt_VSIDS_props_limit*1000000)

    /////////////
    // Statistics

    , dec_vars     (0)
    , decisions    (0)

    ////////////////////
    // Solver references

    , assignmentTrail(s.assignmentTrail)
    , randomNumberGenerator(s.randomNumberGenerator)
    , ca(s.ca)
    , unitPropagator(s.unitPropagator)
    , solver(s)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// STATE MODIFICATION

Lit BranchingHeuristicManager::pickBranchLit() {
    decisions++;

    Var next = var_Undef;
    auto& order_heap = DISTANCE
        ? order_heap_distance
        : ((!VSIDS)? order_heap_CHB:order_heap_VSIDS);

    // Activity based decision:
    while (next == var_Undef || assignmentTrail.value(next) != l_Undef || !decision[next])
        if (order_heap.empty()) {
            return lit_Undef;
        } else {
            handleEventPickBranchLit(solver.conflicts);
            next = order_heap.removeMin();
        }

    return mkLit(next, polarity[next]);
}

void BranchingHeuristicManager::checkSwitchHeuristic(uint64_t propagations) {
    if (switchMode) { 
        switchHeuristic();
        switchMode = false;
    }

    if (propagations - prevPropagations > VSIDS_props_limit){
        prevPropagations = propagations;
        switchMode = true;
        VSIDS_props_limit = VSIDS_props_limit + VSIDS_props_limit / 10;
    }
}

void BranchingHeuristicManager::switchHeuristic(void) {
    VSIDS = !VSIDS;
    if (VSIDS) {
        printf("c Switched to VSIDS.\n");
    } else {
        printf("c Switched to LRB.\n");
    }
    fflush(stdout);

    // Instead of clearing, set vectors to 0
    for (int i = 0; i < assignmentTrail.nVars(); i++) {
        picked[i] = 0;
        conflicted[i] = 0;
        almost_conflicted[i] = 0;
    #ifdef ANTI_EXPLORATION
        canceled[i] = 0;
    #endif
    }
}

static inline void getExponentialDecaySeries(
    vec<int>& decays,
    int length,
    double initialDecay,
    double decayFactor
) {
    while (length--) {
        decays.push(initialDecay);
        initialDecay /= decayFactor;
    }
}

template <class H>
static inline void bumpActivityDistance(
    vec<double>& activity,
    H& heap,
    Var v,
    double bump,
    vec<int>& level_incs
) {
    // Bump variable activity
    activity[v] += bump;

    // Rescale activities if required
    constexpr double RESCALE_THRESHOLD = 1e100;
    if (activity[v] > RESCALE_THRESHOLD) {
        for (Var x = 0; x < activity  .size(); x++) activity  [x] /= RESCALE_THRESHOLD;
        for (int j = 0; j < level_incs.size(); j++) level_incs[j] /= RESCALE_THRESHOLD;
    }
    
    // Update variable's position in priority queue
    if (heap.inHeap(v)) heap.decrease(v);
}

void BranchingHeuristicManager::updateActivityDistance(
    const vec<Var>& involvedVars,
    const vec<int>& varDist,
    int maxDistance
) {
    // Compute exponential decay series
    vec<int> level_incs;
    getExponentialDecaySeries(level_incs, maxDistance, var_iLevel_inc, my_var_decay);

    // Update variable activities
    for (int i = 0; i < involvedVars.size(); i++) {
        Var v = involvedVars[i];
        const double bump = static_cast<double>(varDist[v]) * level_incs[varDist[v] - 1];
        bumpActivityDistance(activity_distance, order_heap_distance, v, bump, level_incs);
    }

    var_iLevel_inc = level_incs.last();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// EVENT HANDLERS

void BranchingHeuristicManager::handleEventLearnedClause(
    const vec<Lit>& learnt_clause,
    vec<bool>& seen,
    int backtrackLevel
) {
    if (VSIDS) {
        // Note: for VSIDS, 'toClear' is set by 'handleEventLitInConflictGraph()'
        for (int i = 0; i < toClear.size(); i++){
            Var v = toClear[i];
            if (assignmentTrail.level(v) >= backtrackLevel - 1)
                varBumpActivity(v, 1);
        }
        toClear.clear();
        varDecayActivity();
    } else {
        // Iterate through every reason clause immediately before the learnt clause
        for (int i = learnt_clause.size() - 1; i >= 0; i--) {
            Var v = var(learnt_clause[i]);
            CRef rea = assignmentTrail.reason(v);
            if (rea == CRef_Undef) continue;
            const Clause& reaC = ca[rea];

            // Iterate through every unique variable in the reason clauses, ignoring variables in
            // the learnt clause
            for (int j = 0; j < reaC.size(); j++){
                Var x = var(reaC[j]);
                if (seen[x]) continue;

                // Increment the 'almost_conflicted' counter
                almost_conflicted[x]++;

                // Mark the variable as seen
                seen[x] = true;
                toClear.push(x);
            }
        }

        // Undo all the changes to seen[]
        for (int j = 0; j < toClear.size(); j++) seen[toClear[j]] = false;
        toClear.clear();
    }
}

void BranchingHeuristicManager::handleEventConflicted(CRef confl, uint64_t conflicts) {
    if (VSIDS) {
        if (--timer == 0 && var_decay < 0.95)
            timer = 5000, var_decay += 0.01;
    } else {
        if (step_size > min_step_size)
            step_size -= step_size_dec;
    }

#if BRANCHING_HEURISTIC == BRANCHING_HEURISTIC_DYNAMIC
    const bool prevDISTANCE = DISTANCE;
    DISTANCE = (conflicts <= 50000);
    if (prevDISTANCE != DISTANCE)
        solver.propagationQueue.prioritizeByActivity(getActivity());
#endif
    if (VSIDS && DISTANCE)
        solver.conflictAnalyzer.collectFirstUIP(confl);
}