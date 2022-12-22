#include "core/BranchingHeuristicManager.h"
#include "core/SolverERTypes.h"
#include "core/Solver.h"

using namespace Minisat;

static const char* _cat = "CORE";

// CHB
static DoubleOption  opt_step_size         (_cat, "step-size",   "Initial step size",                             0.40,     DoubleRange(0, false, 1, false));
static DoubleOption  opt_step_size_dec     (_cat, "step-size-dec","Step size decrement",                          0.000001, DoubleRange(0, false, 1, false));
static DoubleOption  opt_min_step_size     (_cat, "min-step-size","Minimal step size",                            0.06,     DoubleRange(0, false, 1, false));

// VSIDS
static DoubleOption  opt_var_decay         (_cat, "var-decay",   "The variable activity decay factor",            0.80,     DoubleRange(0, false, 1, false));

// Random
static DoubleOption  opt_random_var_freq   (_cat, "rnd-freq",    "The frequency with which the decision heuristic tries to choose a random variable", 0, DoubleRange(0, true, 1, true));
static BoolOption    opt_rnd_init_act      (_cat, "rnd-init",    "Randomize the initial activity", false);

// Phase saving
static IntOption     opt_phase_saving      (_cat, "phase-saving", "Controls the level of phase saving (0=none, 1=limited, 2=full)", 2, IntRange(0, 2));

// Heuristic selection
static IntOption     opt_VSIDS_props_limit (_cat, "VSIDS-lim",  "specifies the number of propagations after which the solver switches between LRB and VSIDS(in millions).", 30, IntRange(1, INT32_MAX));

BranchingHeuristicManager::BranchingHeuristicManager(Solver* s)
#if PRIORITIZE_ER
    : order_heap_VSIDS   (VarOrderLt(activity_VSIDS,    ser->extensionLevel))
    , order_heap_CHB     (VarOrderLt(activity_CHB,      ser->extensionLevel))
    , order_heap_distance(VarOrderLt(activity_distance, ser->extensionLevel))
#else
    : order_heap_VSIDS   (VarOrderLt(activity_VSIDS))
    , order_heap_CHB     (VarOrderLt(activity_CHB))
    , order_heap_distance(VarOrderLt(activity_distance))
#endif

    //////////////////////////
    // Heuristic configuration

    // VSIDS
    , timer    (5000)
    , var_inc  (1)
    , var_decay(opt_var_decay)

    // CHB
    , step_size    (opt_step_size)
    , step_size_dec(opt_step_size_dec)
    , min_step_size(opt_min_step_size)

    // Distance
    , var_iLevel_inc(1)
    , my_var_decay  (0.6)

    // Random
    , random_var_freq(opt_random_var_freq)
    , rnd_pol        (false)
    , rnd_init_act   (opt_rnd_init_act)

    // Phase saving
    , phase_saving(opt_phase_saving)

    // Branching mode
    , VSIDS_props_limit(opt_VSIDS_props_limit*1000000)
    , prev_props (0)
    , switch_mode(false)
    , m_VSIDS    (false)
    , m_DISTANCE (true)

    /////////////
    // Statistics

    , dec_vars     (0)
    , decisions    (0)
    , rnd_decisions(0)

    , conflicts_VSIDS(0)

    ////////////////////
    // Solver references

    , randomNumberGenerator(s->randomNumberGenerator)
    , variableDatabase(s->variableDatabase)
    , ca(s->ca)
    , unitPropagator(s->unitPropagator)
    , solver(s)
{}

BranchingHeuristicManager::~BranchingHeuristicManager() {}

void BranchingHeuristicManager::newVar(Var v, bool sign, bool dvar) {
    // VSIDS
    activity_VSIDS   .push(rnd_init_act ? randomNumberGenerator.drand() * 0.00001 : 0);

    // CHB
    activity_CHB     .push(0);
    conflicted       .push(0);
    almost_conflicted.push(0);
    picked           .push(0);
#ifdef ANTI_EXPLORATION
    canceled         .push(0);
#endif

    // Distance
    activity_distance.push(0);
    pathCs           .push(0);
    var_iLevel       .push(0);
    var_iLevel_tmp   .push(0);

    // Phase saving
    polarity.push(sign);

    // Decision variables
    decision.push();
    setDecisionVar(v, dvar);
}

Lit BranchingHeuristicManager::pickBranchLit() {
    // Update statistics
    decisions++;

    Var next = var_Undef;
    //    Heap<VarOrderLt>& order_heap = VSIDS ? order_heap_VSIDS : order_heap_CHB;
    auto& order_heap = m_DISTANCE ? order_heap_distance : ((!m_VSIDS) ? order_heap_CHB : order_heap_VSIDS);

    // Random decision:
    /*if (randomNumberGenerator.drand() < random_var_freq && !order_heap.empty()){
        next = order_heap[randomNumberGenerator.irand(order_heap.size())];
        if (value(next) == l_Undef && decision[next])
            rnd_decisions++; }*/

    // Activity based decision:
    while (next == var_Undef || variableDatabase.value(next) != l_Undef || !decision[next]) {
        if (order_heap.empty())
            return lit_Undef;
        else {
#ifdef ANTI_EXPLORATION
            if (!m_VSIDS){
                Var v = order_heap_CHB[0];
                uint32_t age = solver->conflicts - canceled[v];
                while (age > 0){
                    double decay = pow(0.95, age);
                    activity_CHB[v] *= decay;
                    if (order_heap_CHB.inHeap(v))
                        order_heap_CHB.increase(v);
                    canceled[v] = solver->conflicts;
                    v = order_heap_CHB[0];
                    age = solver->conflicts - canceled[v];
                }
            }
#endif
            next = order_heap.removeMin();
        }
    }

    return mkLit(next, polarity[next]);
}

void BranchingHeuristicManager::rebuildOrderHeap() {
    vec<Var> vs;
    for (Var v = 0; v < variableDatabase.nVars(); v++)
        if (decision[v] && variableDatabase.value(v) == l_Undef)
            vs.push(v);

    order_heap_VSIDS   .build(vs);
    order_heap_CHB     .build(vs);
    order_heap_distance.build(vs);
}

inline void BranchingHeuristicManager::collectFirstUIPConflictClause(int& minLevel, CRef confl) {
    Clause& c = ca[confl];
    for (int i = 0; i < c.size(); i++) {
        Var v = var(c[i]);
        //        assert(!seen[v]);
        if (solver->level(v) > 0) {
            solver->seen[v] = 1;
            var_iLevel_tmp[v] = 1;
            pathCs[solver->level(v)]++;

            minLevel = std::min(minLevel, solver->level(v));
            assert(minLevel > 0);
        }
    }
}

inline void BranchingHeuristicManager::collectFirstUIPReasonClause(int& minLevel, CRef confl, int reasonVarLevel) {
    Clause& rc = ca[confl];
    // For binary clauses, we don't rearrange literals in propagate(), so check and make sure the first is an implied lit.
    if (rc.size() == 2 && variableDatabase.value(rc[0]) == l_False) {
        assert(variableDatabase.value(rc[1]) != l_False);
        std::swap(rc[0], rc[1]);
    }

    for (int i = 1; i < rc.size(); i++){
        Var v = var(rc[i]);
        assert(solver->level(v) >= 0);

        // Ignore variables set at the root level
        if (solver->level(v) == 0) continue;

        minLevel = std::min(minLevel, solver->level(v));
        if (solver->seen[v]) {
            var_iLevel_tmp[v] = std::min(var_iLevel_tmp[v], static_cast<double>(reasonVarLevel));
        } else {
            var_iLevel_tmp[v] = reasonVarLevel;
            solver->seen[v] = 1;
            pathCs[solver->level(v)]++;
        }
    }
}

// pathCs[k] is the number of variables assigned at level k,
// it is initialized to 0 at the beginning and reset to 0 after the function execution

// Assumes that the trail is in the conflicting state due to confl
void BranchingHeuristicManager::collectFirstUIP(CRef confl) {
    involved_lits.clear();
    
    int minLevel = solver->decisionLevel();
    int maxLevel = 1;

    // Special case for the conflict clause:
    // Increment level path counter for every variable in the clause
    collectFirstUIPConflictClause(minLevel, confl);

    // Iterate backward to at most the first UIP for the conflict clause
    for (int i = solver->trail.size() - 1; i >= solver->trail_lim[minLevel - 1]; i--) {
        Lit p = solver->trail[i]; Var v = var(p);

        // Ignore variables that are not reasons for the conflict
        if (!solver->seen[v]) continue;
        involved_lits.push(p);

        // Since each variable should only appear once on the trail, this unsets
        // the seen flag for this variable for the rest of the function call
        solver->seen[v] = 0;
        
        // Don't explore backward from the variable if it is the first UIP for its decision level
        if (--pathCs[solver->level(v)] == 0) continue;

        // Process reason clause
        const int reasonVarLevel = var_iLevel_tmp[v] + 1;
        collectFirstUIPReasonClause(minLevel, solver->reason(v), reasonVarLevel);
        maxLevel = std::max(maxLevel, reasonVarLevel);
    }

    bumpActivityDistance(involved_lits, maxLevel);
}

void BranchingHeuristicManager::handleEventLearnedClause(const vec<Lit>& out_learnt, const int out_btlevel) {
    if (m_VSIDS) {
        for (int i = 0; i < conflictLits.size(); i++){
            Var v = var(conflictLits[i]);
            if (solver->level(v) >= out_btlevel - 1)
                bumpActivityVSIDS(v, 1);
        }
        conflictLits.clear();
    } else {
        solver->seen[var(out_learnt[0])] = true;
        for (int i = out_learnt.size() - 1; i >= 0; i--) {
            Var v = var(out_learnt[i]);
            CRef rea = solver->reason(v);
            if (rea == CRef_Undef) continue;

            const Clause& reaC = ca[rea];
            for (int i = 0; i < reaC.size(); i++){
                Lit l = reaC[i];
                if (solver->seen[var(l)]) continue;

                solver->seen[var(l)] = true;
                almost_conflicted[var(l)]++;
                solver->analyze_toclear.push(l);
            }
        }
    }

    for (int j = 0; j < solver->analyze_toclear.size(); j++) solver->seen[var(solver->analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)
}

void BranchingHeuristicManager::prioritize(const std::vector<ExtDef>& defs, double scaleFactor) {
    // FIXME: this only forces branching on the last extension variable we add here - maybe add a queue for force branch variables?
    const double desiredActivityCHB   = activity_CHB  [order_heap_CHB  [0]] * scaleFactor;
    const double desiredActivityVSIDS = activity_VSIDS[order_heap_VSIDS[0]] * scaleFactor;
    for (const ExtDef& def : defs) {
        Var v = var(def.x);
        // Prioritize branching on our extension variables
        activity_CHB  [v] = desiredActivityCHB;
        activity_VSIDS[v] = desiredActivityVSIDS;

#if EXTENSION_FORCE_BRANCHING
        // This forces branching because of how branching is implemented when ANTI_EXPLORATION is turned on
        canceled[v] = conflicts;
#endif
        if (order_heap_CHB  .inHeap(v)) order_heap_CHB  .decrease(v);
        if (order_heap_VSIDS.inHeap(v)) order_heap_VSIDS.decrease(v);
    }
}