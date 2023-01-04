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

static const char* _cat = "CORE";

#if BRANCHING_HEURISTIC == VSIDS
// VSIDS
static DoubleOption  opt_var_decay         (_cat, "var-decay",   "The variable activity decay factor",            0.95,     DoubleRange(0, false, 1, false));
#endif

#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
// CHB
static DoubleOption  opt_step_size         (_cat, "step-size",     "Initial step size",   0.40,     DoubleRange(0, false, 1, false));
static DoubleOption  opt_step_size_dec     (_cat, "step-size-dec", "Step size decrement", 0.000001, DoubleRange(0, false, 1, false));
static DoubleOption  opt_min_step_size     (_cat, "min-step-size", "Minimal step size",   0.06,     DoubleRange(0, false, 1, false));
#endif
#if BRANCHING_HEURISTIC == CHB
static DoubleOption  opt_reward_multiplier (_cat, "reward-multiplier", "Reward multiplier", 0.9, DoubleRange(0, true, 1, true));
#endif

// Random
static DoubleOption  opt_random_var_freq   (_cat, "rnd-freq", "The frequency with which the decision heuristic tries to choose a random variable", 0, DoubleRange(0, true, 1, true));
static BoolOption    opt_rnd_init_act      (_cat, "rnd-init", "Randomize the initial activity", false);

// Phase saving
static IntOption     opt_phase_saving      (_cat, "phase-saving", "Controls the level of phase saving (0=none, 1=limited, 2=full)", 2, IntRange(0, 2));

// Heuristic selection
static IntOption     opt_VSIDS_props_limit (_cat, "VSIDS-lim", "specifies the number of propagations after which the solver switches between LRB and VSIDS(in millions).", 30, IntRange(1, INT32_MAX));

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

BranchingHeuristicManager::BranchingHeuristicManager(Solver& s)
#if PRIORITIZE_ER
#ifdef EXTLVL_ACTIVITY
    : order_heap(
        VarOrderLt<double>{extensionLevelActivity, false},
        VarOrderLt<double>{activity, false},
        extensionLevel
    )
#else
    : order_heap_extlvl(VarOrderLt2<uint64_t, double>{extensionLevel, true, activity, false})
    , order_heap_degree(VarOrderLt2<uint64_t, double>{degree, true, activity, false})
#endif
#else
    : order_heap(VarOrderLt<double>{activity, false})
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

#if BRANCHING_HEURISTIC == CHB
    , prev_trail_index(0)
    , reward_multiplier(opt_reward_multiplier)
#endif

    // Random
    , random_var_freq(opt_random_var_freq)
    , rnd_pol        (false)
    , rnd_init_act   (opt_rnd_init_act)

    // Phase saving
    , phase_saving(static_cast<PhaseSavingLevel>(static_cast<int>(opt_phase_saving)))

    /////////////
    // Statistics

    , dec_vars     (0)
    , decisions    (0)
    , rnd_decisions(0)

    ////////////////////
    // Solver references

    , assignmentTrail(s.assignmentTrail)
    , randomNumberGenerator(s.randomNumberGenerator)
    , ca(s.ca)
    , unitPropagator(s.unitPropagator)
    , solver(s)
{
#ifdef POLARITY_VOTING
    group_polarity.push(0);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// STATE MODIFICATION

Lit BranchingHeuristicManager::pickBranchLit() {
    decisions++;
    Var next = var_Undef;

    // Random decision:
    {
    #if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
        auto& order_heap = order_heap_extlvl.empty() ? order_heap_degree : order_heap_extlvl;
    #endif
        if (randomNumberGenerator.drand() < random_var_freq && !order_heap.empty()){
            next = order_heap[randomNumberGenerator.irand(order_heap.size())];
            if (assignmentTrail.value(next) == l_Undef && decision[next])
                rnd_decisions++;
        }
    }

    // Activity based decision:
    while (next == var_Undef || assignmentTrail.value(next) != l_Undef || !decision[next]) {
    #if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
        auto& order_heap = order_heap_extlvl.empty() ? order_heap_degree : order_heap_extlvl;
    #endif
        if (order_heap.empty()){
            next = var_Undef;
            break;
        } else {
        #if BRANCHING_HEURISTIC == LRB
            handleEventPickBranchLitLRB(solver.conflicts);
        #endif
            next = order_heap.removeMin();
        }
    }

    // No literals remaining! Skip polarity selection
    if (next == var_Undef) return lit_Undef;

    // Check option to select a random polarity
    if (rnd_pol) return mkLit(next, randomNumberGenerator.drand() < 0.5);

#ifdef POLARITY_VOTING
    // Vote for the next branch literal
    const double vote = group_polarity[extensionLevel[next]];
    const bool preferred_polarity = (vote == 0) ? polarity[next] : (vote < 0);

    // Update stats
    group_polarity[extensionLevel[next]] += (preferred_polarity ? (-1) : (+1)) * 0.01;
#else
    const bool preferred_polarity = polarity[next];
#endif

    // Return a literal from the selected variable and polarity
    return mkLit(next, preferred_polarity);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// VARIABLE SELECTION HEURISTIC STATE MODIFICATION

void BranchingHeuristicManager::rebuildPriorityQueue() {
    const int numVars = assignmentTrail.nVars();

    // Build the priority queue
#if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
    order_heap_extlvl.clear();
    order_heap_degree.clear();
    for (Var v = 0; v < numVars; v++) {
        if (decision[v] && assignmentTrail.value(v) == l_Undef) {
            order_heap_degree.insert(v);
            if (extensionLevel[v]) order_heap_extlvl.insert(v);
        }
    }
#else
    vec<Var> vs;
    for (Var v = 0; v < numVars; v++) {
        if (decision[v] && assignmentTrail.value(v) == l_Undef) {
            vs.push(v);
        }
    }
    order_heap.build(vs);
#endif
}

void BranchingHeuristicManager::handleEventLearnedClause(
    const vec<Lit>& learnt_clause,
    vec<bool>& seen
) {
#if ALMOST_CONFLICT
    // Iterate through every reason clause immediately before the learnt clause
    for (int i = learnt_clause.size() - 1; i >= 0; i--) {
        Var v = var(learnt_clause[i]);
        CRef rea = assignmentTrail.reason(v);
        if (rea == CRef_Undef) continue;
        Clause& reaC = ca[rea];

        // Iterate through every unique variable in the reason clauses, ignoring
        // variables in the learnt clause
        for (int j = 0; j < reaC.size(); j++) {
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
#endif

#ifdef POLARITY_VOTING
    // Apply EMA to group polarities
    for (int k = 0; k < group_polarity.size(); k++) {
        if (polarity_count[k]) {
            group_polarity[k] = 0.9 * (group_polarity[k] + polarity_count[k]);
        }
    }
#endif

#if PRIORITIZE_ER
    for (int k = 0; k < learnt_clause.size(); k++)
        degree[var(learnt_clause[k])]++;
#endif
}