#include "core/BranchingHeuristicManager.h"
#include "core/Solver.h"

using namespace Minisat;

static const char* _cat = "CORE";

#if BRANCHING_HEURISTIC == VSIDS
// VSIDS
static DoubleOption  opt_var_decay         (_cat, "var-decay",   "The variable activity decay factor",            0.95,     DoubleRange(0, false, 1, false));
#endif

#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
// CHB
static DoubleOption  opt_step_size         (_cat, "step-size",   "Initial step size",                             0.40,     DoubleRange(0, false, 1, false));
static DoubleOption  opt_step_size_dec     (_cat, "step-size-dec","Step size decrement",                          0.000001, DoubleRange(0, false, 1, false));
static DoubleOption  opt_min_step_size     (_cat, "min-step-size","Minimal step size",                            0.06,     DoubleRange(0, false, 1, false));
#endif

// Random
static DoubleOption  opt_random_var_freq   (_cat, "rnd-freq",    "The frequency with which the decision heuristic tries to choose a random variable", 0, DoubleRange(0, true, 1, true));
static BoolOption    opt_rnd_init_act      (_cat, "rnd-init",    "Randomize the initial activity", false);

// Phase saving
static IntOption     opt_phase_saving      (_cat, "phase-saving", "Controls the level of phase saving (0=none, 1=limited, 2=full)", 2, IntRange(0, 2));

// Heuristic selection
static IntOption     opt_VSIDS_props_limit (_cat, "VSIDS-lim",  "specifies the number of propagations after which the solver switches between LRB and VSIDS(in millions).", 30, IntRange(1, INT32_MAX));

BranchingHeuristicManager::BranchingHeuristicManager(Solver* s)
#if PRIORITIZE_ER
#ifdef EXTLVL_ACTIVITY
  : order_heap         (VarOrderLt(extensionLevelActivity), VarOrderLt(activity), extensionLevel)
#else
  : order_heap_extlvl  (VarOrderLt(activity, extensionLevel, false))
  , order_heap_degree  (VarOrderLt(activity, degree, true))
#endif
#else
  : order_heap         (VarOrderLt(activity))
#endif

    //////////////////////////
    // Heuristic configuration

    // VSIDS
#if BRANCHING_HEURISTIC == VSIDS
    , var_inc            (1)
    , var_decay(opt_var_decay)
#endif

    // CHB
#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
    , step_size    (opt_step_size)
    , step_size_dec(opt_step_size_dec)
    , min_step_size(opt_min_step_size)
#endif

    // Random
    , random_var_freq(opt_random_var_freq)
    , rnd_pol        (false)
    , rnd_init_act   (opt_rnd_init_act)

    // Phase saving
    , phase_saving(opt_phase_saving)

    /////////////
    // Statistics

    , dec_vars     (0)
    , decisions    (0)
    , rnd_decisions(0)

    ////////////////////
    // Solver references

    , assignmentTrail(s->assignmentTrail)
    , randomNumberGenerator(s->randomNumberGenerator)
    , variableDatabase(s->variableDatabase)
    , ca(s->ca)
    , unitPropagator(s->unitPropagator)
    , solver(s)
{
#ifdef POLARITY_VOTING
    group_polarity.push(0);
#endif
}

BranchingHeuristicManager::~BranchingHeuristicManager() {}

void BranchingHeuristicManager::newVar(Var v, bool sign, bool dvar) {
    // VSIDS
    activity.push(rnd_init_act ? randomNumberGenerator.drand() * 0.00001 : 0);

    // CHB
    conflicted.push(0);
#if ALMOST_CONFLICT
    almost_conflicted.push(0);
#endif
    picked.push(0);
#if ANTI_EXPLORATION
    canceled.push(0);
#endif

    // Phase saving
    polarity.push(sign);

    // Decision variables
    decision.push();
    setDecisionVar(v, dvar);

#if PRIORITIZE_ER
    degree.push(0);
    extensionLevel.push(0);
#endif
#if BRANCHING_HEURISTIC == CHB
    last_conflict.push(0);
#endif

    // Statistics
    total_actual_rewards.push(0);
    total_actual_count.push(0);
}

Lit BranchingHeuristicManager::pickBranchLit() {
    decisions++;
    Var next = var_Undef;

    {
#if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
    Heap<VarOrderLt>& order_heap = order_heap_extlvl.empty() ? order_heap_degree : order_heap_extlvl;
#endif

    // Random decision:
    if (randomNumberGenerator.drand() < random_var_freq && !order_heap.empty()){
        next = order_heap[randomNumberGenerator.irand(order_heap.size())];
        if (variableDatabase.value(next) == l_Undef && decision[next])
            rnd_decisions++; }
    }

    // Activity based decision:
    while (next == var_Undef || variableDatabase.value(next) != l_Undef || !decision[next]) {
#if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
        Heap<VarOrderLt>& order_heap = order_heap_extlvl.empty() ? order_heap_degree : order_heap_extlvl;
#endif
        if (order_heap.empty()){
            next = var_Undef;
            break;
        } else {
#if ANTI_EXPLORATION
            next = order_heap[0];
            uint64_t age = solver->conflicts - canceled[next];
            while (age > 0 && variableDatabase.value(next) == l_Undef) {
                double decay = pow(0.95, age);
                activity[next] *= decay;
                if (order_heap.inHeap(next)) {
                    order_heap.increase(next);
                }
                canceled[next] = solver->conflicts;
                next = order_heap[0];
                age = solver->conflicts - canceled[next];
            }
#endif
            next = order_heap.removeMin();
        }
    }

    // No literals remaining! Skip polarity selection
    if (next == var_Undef) return lit_Undef;

    if (rnd_pol) return mkLit(next, randomNumberGenerator.drand() < 0.5);

#ifdef POLARITY_VOTING
    // Vote for the next branch literal
    double vote = group_polarity[extensionLevel[next]];
    bool preferred_polarity = (vote == 0) ? polarity[next] : (vote < 0);
#else
    bool preferred_polarity = polarity[next];
#endif

#ifdef POLARITY_VOTING
    // Update stats
    group_polarity[extensionLevel[next]] += (preferred_polarity ? (-1) : (+1)) * 0.01;
#endif
    return mkLit(next, preferred_polarity);
}


void BranchingHeuristicManager::rebuildOrderHeap() {
#if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
    order_heap_extlvl.clear();
    order_heap_degree.clear();
    for (Var v = 0; v < variableDatabase.nVars(); v++)
        if (decision[v] && variableDatabase.(v) == l_Undef) {
            order_heap_degree.insert(v);
            if (extensionLevel[v]) order_heap_extlvl.insert(v);
        }
#else
    vec<Var> vs;
    for (Var v = 0; v < variableDatabase.nVars(); v++)
        if (decision[v] && variableDatabase.value(v) == l_Undef)
            vs.push(v);
    order_heap.build(vs);
#endif
}

void BranchingHeuristicManager::handleEventLearnedClause(const vec<Lit>& out_learnt, const int out_btlevel) {
#if ALMOST_CONFLICT
    for(int i = out_learnt.size() - 1; i >= 0; i--) {
        Var v = var(out_learnt[i]);
        CRef rea = assignmentTrail.reason(v);
        if (rea != CRef_Undef) {
            Clause& reaC = ca[rea];
            for (int i = 0; i < reaC.size(); i++) {
                Lit l = reaC[i];
                if (!solver->seen[var(l)]) {
                    solver->seen[var(l)] = true;
                    almost_conflicted[var(l)]++;
                    solver->analyze_toclear.push(l);
                }
            }
        }
    }
#endif
    for (int j = 0; j < solver->analyze_toclear.size(); j++) solver->seen[var(solver->analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)

#ifdef POLARITY_VOTING
    // Apply EMA to group polarities
    for (int k = 0; k < group_polarity.size(); k++) {
        if (polarity_count[k]) {
            group_polarity[k] = 0.9 * (group_polarity[k] + polarity_count[k]);
        }
    }
#endif
}