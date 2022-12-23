/****************************************************************************************[Solver.h]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010, Niklas Sorensson

Chanseok Oh's MiniSat Patch Series -- Copyright (c) 2015, Chanseok Oh
 
Maple_LCM, Based on MapleCOMSPS_DRUP -- Copyright (c) 2017, Mao Luo, Chu-Min LI, Fan Xiao: implementing a learnt clause minimisation approach
Reference: M. Luo, C.-M. Li, F. Xiao, F. Manya, and Z. L. , “An effective learnt clause minimization approach for cdcl sat solvers,” in IJCAI-2017, 2017, pp. to–appear.

Maple_LCM_Dist, Based on Maple_LCM -- Copyright (c) 2017, Fan Xiao, Chu-Min LI, Mao Luo: using a new branching heuristic called Distance at the beginning of search 

xMaple_LCM_Dist, based on Maple_LCM_Dist -- Copyright (c) 2022, Jonathan Chung, Vijay Ganesh, Sam Buss

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

#ifndef Minisat_BranchingHeuristicManager_h
#define Minisat_BranchingHeuristicManager_h

#include <math.h>
#include "core/SolverTypes.h"
#include "core/RandomNumberGenerator.h"
#include "core/VariableDatabase.h"
#include "core/UnitPropagator.h"
#include "mtl/Heap.h"

namespace Minisat {
    // Forward declarations
    class Solver;

    /**
     * @brief This class handles variable branching.
     * 
     */
    class BranchingHeuristicManager {
    public:
    protected:
        // Comparator for priority queue
        struct VarOrderLt {
            const vec<double>& activity;
#if PRIORITIZE_ER
            const vec<unsigned int>& extensionLevel;
            bool operator () (Var x, Var y) const {
#if PRIORITIZE_ER_LOW
                if (extensionLevel[x] != extensionLevel[y]) return extensionLevel[x] < extensionLevel[y];
#else
                if (extensionLevel[x] != extensionLevel[y]) return extensionLevel[x] > extensionLevel[y];
#endif
                else                                        return activity[x] > activity[y];
            }
            VarOrderLt(const vec<double>& act, const vec<unsigned int>& extlvl)
                : activity(act)
                , extensionLevel(extlvl)
            { }
#else
            bool operator () (Var x, Var y) const { return activity[x] > activity[y]; }
            VarOrderLt(const vec<double>&  act) : activity(act) { }
#endif
        };

        //////////////////////
        // MEMBER VARIABLES //
        //////////////////////

        // A heuristic measurement of the activity of a variable.
        vec<double> activity;

    #if PRIORITIZE_ER
    #ifdef EXTLVL_ACTIVITY
        Multiheap<VarOrderLt> order_heap;
    #else
        Heap<VarOrderLt>    order_heap_extlvl;       // A priority queue of variables ordered with respect to the extension level.
        Heap<VarOrderLt>    order_heap_degree;       // A priority queue of variables ordered with respect to the variable degree.
    #endif
    #else
        Heap<VarOrderLt>    order_heap;       // A priority queue of variables ordered with respect to the variable activity.
    #endif

        vec<char> decision; // Declares whether a variable is eligible for selection in the decision heuristic.
        
        //////////////////////////
        // Heuristic configuration

        // VSIDS
#if BRANCHING_HEURISTIC == VSIDS
        double var_inc; // Amount to bump next variable with.
        double var_decay;
        vec<Lit> conflictLits; // Literals that participate in the conflict graph
#endif

        // CHB
#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
        double step_size;
        double step_size_dec;
        double min_step_size;
#endif
        vec<uint64_t> conflicted;
#if ALMOST_CONFLICT
        vec<uint64_t> almost_conflicted;
#endif
        vec<uint64_t> picked;
#ifdef ANTI_EXPLORATION
        vec<uint64_t> canceled;
#endif

#if BRANCHING_HEURISTIC == CHB
        vec<uint64_t> last_conflict;
        int action;
        double reward_multiplier;
#endif

        // Random
        double random_var_freq;
        bool   rnd_pol;      // Use random polarities for branching heuristics.
        bool   rnd_init_act; // Initialize variable activities with a small random value.

        // Phase saving
        int phase_saving;   // Controls the level of phase saving (0=none, 1=limited, 2=full).
        vec<char> polarity; // The preferred polarity of each variable.

    #if PRIORITIZE_ER
        // Number of times a variable occurs in a clause
        vec<uint64_t> degree;
        // Map from variables to their extension level
        vec<uint64_t> extensionLevel;

    #ifdef EXTLVL_ACTIVITY
        vec<double> extensionLevelActivity; // The activity of each extension level
    #endif
    #ifdef POLARITY_VOTING
        vec<unsigned int> polarity_count;
        vec<double> group_polarity;   // The preferred polarity of each group.
    #endif
    #endif

    public:
        ////////////////
        // STATISTICS //
        ////////////////

        uint64_t dec_vars     ; // Total number of decision variables
        uint64_t decisions    ; // Total number of decisions
        uint64_t rnd_decisions; // Total number of random decisions

        vec<long double> total_actual_rewards;
        vec<int> total_actual_count;

    protected:
        ///////////////////////
        // SOLVER REFERENCES //
        ///////////////////////

        AssignmentTrail& assignmentTrail;
        RandomNumberGenerator& randomNumberGenerator;
        VariableDatabase& variableDatabase;
        ClauseAllocator& ca;
        UnitPropagator& unitPropagator;
        Solver* solver;

    public:

        /**
         * @brief Construct a new BranchingHeuristicManager object
         * 
         * @param s Pointer to main solver object - must not be nullptr
         */
        BranchingHeuristicManager(Solver* s);

        /**
         * @brief Destroy the BranchingHeuristicManager object
         * 
         */
        ~BranchingHeuristicManager();

        /**
         * @brief Set up internal data structures for a new variable
         * 
         * @param v the variable to register
         * @note Assumes that @code{v} has not already been registered
         */
        void newVar(Var v, bool sign, bool dvar);

        /**
         * @brief Insert a variable in the decision order priority queue.
         * 
         * @param x the variable to insert
         */
        void insertVarOrder(Var x);

        /**
         * @brief Return the next decision variable.
         * 
         * @return the decision literal
         */
        Lit pickBranchLit();

        /**
         * @brief Rebuild the priority queue from scratch
         * 
         */
        void rebuildOrderHeap();

        ////////////////
        // HEURISTICS //
        ////////////////

#if BRANCHING_HEURISTIC == VSIDS
        // VSIDS

        void decayActivityVSIDS();                   // Decay all variables with the specified factor. Implemented by increasing the 'bump' value instead.
        void bumpActivityVSIDS (Var v, double mult); // Increase a variable with the current 'bump' value.
#endif

        // Phase saving

        void setPolarity   (Var v, bool b); // Declare which polarity the decision heuristic should use for a variable. Requires mode 'polarity_user'.
        void setDecisionVar(Var v, bool b); // Declare if a variable should be eligible for selection in the decision heuristic.

        /////////////////////
        // EVENT LISTENERS //
        /////////////////////

        /**
         * @brief Update data structures for branching heuristics upon variable assignment.
         * 
         * @param l The literal that was assigned
         * @param conflicts The current number of conflicts in the solver
         */
        void handleEventLitAssigned(Lit l, uint64_t conflicts);
        
        /**
         * @brief Update data structures for branching heuristics upon variable unassignment.
         * 
         * @param l The literal that was unassigned
         * @param conflicts The current number of conflicts in the solver
         * @param assignedAtLastLevel true if l was assigned at the last decision level, false otherwise
         */
        void handleEventLitUnassigned(Lit l, uint64_t conflicts, bool assignedAtLastLevel);

        /**
         * @brief Update data structures for branching heuristics when unit propagation exits.
         * 
         * @param conflicts the new total number of conflicts
         * @param conflicted true iff propagation resulted in conflict
         */
        void handleEventPropagated(uint64_t conflicts, bool conflicted);

        /**
         * @brief Update data structures for branching heuristics upon conflict.
         * 
         * @param conflicts the new total number of conflicts
         */
        void handleEventConflicted(uint64_t conflicts);

        /**
         * @brief Update data structures for branching heuristics when a literal appears in
         * the conflict graph.
         * 
         * @param l The literal in the conflict graph
         * @param conflicts the previous total number of conflicts
         */
        void handleEventLitInConflictGraph(Lit l, uint64_t conflicts);

        /**
         * @brief Update data structures for branching heuristics when the solver restarts.
         * 
         * @param l The literal in the conflict graph
         */
        void handleEventRestarted(uint64_t propagations);

        /**
         * @brief Update data structures for branching heuristics after learning a clause
         * 
         * @param out_learnt  the learnt clause
         * @param out_btlevel the backtrack level for the learnt clause
         */
        void handleEventLearnedClause(const vec<Lit>& out_learnt, const int out_btlevel);

        ///////////////
        // ACCESSORS //
        ///////////////

        const vec<double>& getActivityVSIDS() const;

    private:
        /////////////////////////////////////
        // HELPER FUNCTIONS FOR HEURISTICS //
        /////////////////////////////////////

    };

    inline void BranchingHeuristicManager::insertVarOrder(Var x) {
    #if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
        Heap<VarOrderLt>& order_heap = extensionLevel[x] ? order_heap_extlvl : order_heap_degree;
    #endif
        if (!order_heap.inHeap(x) && decision[x]) order_heap.insert(x);
    }

#if BRANCHING_HEURISTIC == VSIDS
    inline void BranchingHeuristicManager::decayActivityVSIDS() { var_inc *= (1 / var_decay); }

    inline void BranchingHeuristicManager::bumpActivityVSIDS(Var v, double mult) {
        if ((activity[v] += var_inc) > 1e100) {
    #if PRIORITIZE_ER && defined(EXTLVL_ACTIVITY)
            // Clear extension level activity
            for (int i = 0; i < extensionLevelActivity.size(); i++) {
                extensionLevelActivity[i] = 0;
            }
    #endif
            // Rescale:
            for (int i = 0; i < variableDatabase.nVars(); i++) {
                activity[i] *= 1e-100;
    #if PRIORITIZE_ER && defined(EXTLVL_ACTIVITY)
                extensionLevelActivity[extensionLevel[i]] += activity[i];
    #endif
            }
            var_inc *= 1e-100;
        }

        // Update order_heap with respect to new activity:
    #if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
        if (order_heap_extlvl.inHeap(v))
            order_heap_extlvl.decrease(v);
        if (order_heap_degree.inHeap(v))
            order_heap_degree.decrease(v);
    #else
        if (order_heap.inHeap(v))
            order_heap.decrease(v);
    #endif
    }
#endif

    inline void BranchingHeuristicManager::setPolarity   (Var v, bool b) { polarity[v] = b; }
    inline void BranchingHeuristicManager::setDecisionVar(Var v, bool b) { 
        if      ( b && !decision[v]) dec_vars++;
        else if (!b &&  decision[v]) dec_vars--;

        decision[v] = b;
        insertVarOrder(v);
    }

    inline void BranchingHeuristicManager::handleEventLitAssigned(Lit l, uint64_t conflicts) {
        const Var v = var(l);
        assert(variableDatabase.value(l) == l_Undef);
        picked[v] = conflicts;
    #if ANTI_EXPLORATION
        uint64_t age = conflicts - canceled[v];
        if (age > 0) {
            double decay = pow(0.95, age);
            activity[v] *= decay;
    #if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
            if (order_heap_extlvl.inHeap(v))
                order_heap_extlvl.increase(v);
            if (order_heap_degree.inHeap(v))
                order_heap_degree.increase(v);
    #else
            if (order_heap.inHeap(v)) {
                order_heap.increase(v);
            }
    #endif
        }
    #endif
        conflicted[v] = 0;
    #if ALMOST_CONFLICT
        almost_conflicted[v] = 0;
    #endif
    }

    inline void BranchingHeuristicManager::handleEventLitUnassigned(Lit l, uint64_t conflicts, bool assignedAtLastLevel) {
        const Var v = var(l);

        uint64_t age = conflicts - picked[v];
        if (age > 0) {
            double reward = ((double) conflicted[v]) / ((double) age);
#if BRANCHING_HEURISTIC == LRB
#if ALMOST_CONFLICT
            double adjusted_reward = ((double) (conflicted[v] + almost_conflicted[v])) / ((double) age);
#else
            double adjusted_reward = reward;
#endif
            double old_activity = activity[v];
            activity[v] = step_size * adjusted_reward + ((1 - step_size) * old_activity);
#if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
            auto& order_heap = solver->extensionLevel[v] ? solver->order_heap_extlvl : solver->order_heap_degree;
#endif
            if (order_heap.inHeap(v)) {
                if (activity[v] > old_activity)
                    order_heap.decrease(v);
                else
                    order_heap.increase(v);
            }
#endif

            // Update statistics
            total_actual_rewards[v] += reward;
            total_actual_count[v] ++;
        }
#if ANTI_EXPLORATION
        canceled[v] = conflicts;
#endif

        // Phase saving
        if (phase_saving > 1 || (phase_saving == 1) && assignedAtLastLevel)
            setPolarity(v, sign(l));

        // Update priority queue
        insertVarOrder(v);
    }

    inline void BranchingHeuristicManager::handleEventPropagated(uint64_t conflicts, bool conflicted) {
    #if BRANCHING_HEURISTIC == CHB
        const double multiplier = conflicted ? reward_multiplier : 1.0;
        for (int a = action; a < assignmentTrail.nAssigns(); a++) {
            Var v = var(assignmentTrail[a]);
            const uint64_t age = conflicts - last_conflict[v] + 1;
            const double reward = multiplier / age ;
            const double old_activity = activity[v];
            activity[v] = step_size * reward + ((1 - step_size) * old_activity);
            if (order_heap.inHeap(v)) {
                if (activity[v] > old_activity)
                    order_heap.decrease(v);
                else
                    order_heap.increase(v);
            }
        }

        action = assignmentTrail.nAssigns();
    #endif
    }

    inline void BranchingHeuristicManager::handleEventConflicted(uint64_t conflicts) {
    #if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
        if (step_size > min_step_size)
            step_size -= step_size_dec;
    #endif

    #ifdef POLARITY_VOTING
        // Count votes for the polarity that led to the conflict
        polarity_count.clear();
        for (int k = 0; k < group_polarity.size(); k++) polarity_count.push(0);
    #endif
    }

    inline void BranchingHeuristicManager::handleEventLitInConflictGraph(Lit q, uint64_t conflicts) {
    #if BRANCHING_HEURISTIC == CHB
        last_conflict[var(q)] = conflicts;
    #elif BRANCHING_HEURISTIC == VSIDS
        bumpActivityVSIDS(var(q), var_inc);
    #endif
    #ifdef POLARITY_VOTING
        // Count votes for the polarity that led to the conflict
        count[extensionLevel[var(q)]] += sign(q) ? (+1) : (-1);
    #endif
        conflicted[var(q)]++;
    }

    inline void BranchingHeuristicManager::handleEventRestarted(uint64_t propagations) {

    }

    inline const vec<double>& BranchingHeuristicManager::getActivityVSIDS() const { return activity; }
} // namespace Minisat

#endif