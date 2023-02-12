/*********************************************************************[BranchingHeuristicManager.h]
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

#ifndef Minisat_BranchingHeuristicManager_h
#define Minisat_BranchingHeuristicManager_h

#include <math.h>
#include "core/AssignmentTrail.h"
#include "core/RandomNumberGenerator.h"
#include "mtl/Heap.h"

namespace Minisat {
    // Forward declarations
    class Solver;
    class UnitPropagator;

    /**
     * @brief This class selects decision variables for branching and manages data structures for
     * branching heuristics.
     * 
     */
    class BranchingHeuristicManager {
    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // HELPER TYPES

        /**
         * @brief Comparator for priority queue - orders variables according to an activity metric
         * 
         */
        template<typename T>
        struct VarOrderLt {
            /// @brief The metric for the sorting order
            const vec<T>& activity;
            bool ascending;

            /**
             * @brief Comparison operator: compare two variables and determine whether one should
             * be sorted before the other
             * 
             * @param x The first variable to compare
             * @param y The second variable to compare
             * @return true iff x should be sorted before y
             */
            bool operator () (Var x, Var y) const {
                return ascending ^ (activity[x] > activity[y]);
            }
        };

        /**
         * @brief Comparator for priority queue - orders variables according to two activity
         * metrics
         * 
         */
        template<typename T1, typename T2>
        struct VarOrderLt2 {
            /// @brief The metric for the greatest-priority sorting order
            const vec<T1>& activity1;

            /// @brief True iff @code{activity1} should be sorted in ascending order 
            const bool ascending1;

            /// @brief The metric for the least-priority sorting order
            const vec<T2>& activity2;

            /// @brief True iff @code{activity2} should be sorted in ascending order
            const bool ascending2;

            /**
             * @brief Comparison operator: compare two variables and determine whether one should
             * be sorted before the other
             * 
             * @param x The first variable to compare
             * @param y The second variable to compare
             * @return true iff x should be sorted before y
             */
            bool operator () (Var x, Var y) const {
                if (activity1[x] != activity1[y])
                    return ascending1 ^ (activity1[x] > activity1[y]);
                else
                    return ascending2 ^ (activity2[x] > activity2[y]);
            }
        };

        enum PhaseSavingLevel: int {
            NONE    = 0,
            LIMITED = 1,
            FULL    = 2,
        };

    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // VARIABLE SELECTION HEURISTIC MEMBER VARIABLES

        /// @brief Declares whether a variable is eligible for selection in the decision heuristic.
        vec<bool> decision;

        /// @brief A heuristic measurement of the activity of a variable.
        vec<double> activity;

    #if PRIORITIZE_ER
        /// @brief Number of times a variable occurs in a clause
        vec<uint64_t> degree;

        /// @brief The extension level of each variable
        vec<uint64_t> extensionLevel;

    #ifdef EXTLVL_ACTIVITY
        /// @brief The activity of each extension level
        vec<double> extensionLevelActivity;

        /// @brief A priority queue of variables ordered with respect to the extension level
        /// activity.
        Multiheap<VarOrderLt> order_heap;
    #else
        /// @brief A priority queue of variables ordered with respect to the extension level.
        Heap< VarOrderLt2<uint64_t, double> > order_heap_extlvl;
        
        /// @brief A priority queue of variables ordered with respect to the variable degree.
        Heap< VarOrderLt2<uint64_t, double> > order_heap_degree;
    #endif
    #else
        /// @brief A priority queue of variables ordered with respect to the variable activity.
        Heap< VarOrderLt<double> > order_heap;
    #endif

    #if BRANCHING_HEURISTIC == CHB
        ////////
        // CHB

        /// @brief The number of times a variable appears in a conflict graph
        vec<uint64_t> conflicted;

        /// @brief The number of conflicts seen by the solver before a variable was assigned
        vec<uint64_t> picked;

        /// @brief The number of conflicts seen by the solver before a variable occurred in a
        /// conflict graph
        vec<uint64_t> last_conflict;

    #elif BRANCHING_HEURISTIC == LRB
        ////////
        // LRB

        /// @brief The number of times a variable appears in a conflict graph
        vec<uint64_t> conflicted;

        /// @brief The total number of conflicts seen by the solver before the variable was
        /// assigned
        vec<uint64_t> picked;

    #if ALMOST_CONFLICT
        /// @brief The number of times a variable appeared directly before the cut in the conflict
        /// graph for a learnt clause
        vec<uint64_t> almost_conflicted;
    #endif
    #if ANTI_EXPLORATION
        /// @brief The total number of conflicts seen by the solver before the variable was
        /// selected as a decision literal or unassigned
        vec<uint64_t> canceled;
    #endif
    #endif

        ///////////////////////////////////////////////////////////////////////////////////////////
        // POLARITY SELECTION HEURISTIC MEMBER VARIABLES

        /// @brief The preferred polarity of each variable.
        vec<char> polarity;

    #if PRIORITIZE_ER && defined(POLARITY_VOTING)
        vec<unsigned int> polarity_count;
        
        // The preferred polarity of each group.
        vec<double> group_polarity;
    #endif
        
    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // PARAMETERS

    #if BRANCHING_HEURISTIC == VSIDS
        //////////
        // VSIDS

        /// @brief Amount by which to bump variable activity
        double var_inc;

        /// @brief Amount by which to decay variable activities
        double var_decay;

    #elif BRANCHING_HEURISTIC == CHB
        ////////
        // CHB

        /// @brief Step size for computing CHB activity
        double step_size;

        /// @brief Step size decrement
        double step_size_dec;

        /// @brief Minimum step size
        double min_step_size;

        /// @brief Total number of variable assignments after the previous call to BCP
        int prev_trail_index;

        /// @brief CHB reward multiplier
        double reward_multiplier;

    #elif BRANCHING_HEURISTIC == LRB
        ////////
        // LRB

        /// @brief Step size for computing LRB activity
        double step_size;

        /// @brief Step size decrement
        double step_size_dec;

        /// @brief Minimum step size
        double min_step_size;
    #endif

        ///////////
        // Random

        /// @brief The probabilistic frequency for selecting variables at random
        double random_var_freq;

        /// @brief Use random polarities for branching heuristics
        bool rnd_pol;

        /// @brief Initialize variable activities with a small random value
        bool rnd_init_act;

        /////////////////
        // Phase saving

        /// @brief Controls the level of phase saving (0=none, 1=limited, 2=full).
        PhaseSavingLevel phase_saving;

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // TEMPORARY VARIABLES

        /// @brief Temporary list of variables whose values need to be cleared
        vec<Var> toClear;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATISTICS

        /// @brief Total number of decision variables
        uint64_t dec_vars;

        /// @brief Total number of decisions
        uint64_t decisions;

        /// @brief Total number of random decisions
        uint64_t rnd_decisions;

    #if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
        /// @brief Total lifetime rewards for each variable
        vec<long double> total_actual_rewards;

        /// @brief Total number of times the variable received a reward
        vec<int> total_actual_count;
    #endif

    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // SOLVER REFERENCES

        AssignmentTrail& assignmentTrail;
        RandomNumberGenerator& randomNumberGenerator;
        ClauseAllocator& ca;
        UnitPropagator& unitPropagator;
        Solver& solver;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CONSTRUCTORS

        /**
         * @brief Construct a new BranchingHeuristicManager object
         * 
         * @param s Reference to main solver object
         */
        BranchingHeuristicManager(Solver& s);

        /**
         * @brief Destroy the BranchingHeuristicManager object
         * 
         */
        ~BranchingHeuristicManager() = default;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // ACCESSORS

        /**
         * @brief Get the VSIDS activity array
         * 
         * @return the VSIDS activity array 
         */
        const vec<double>& getActivityVSIDS(void) const;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATE MODIFICATION

        /**
         * @brief Set up internal data structures for a new variable
         * 
         * @param v the variable to register
         * @note Assumes that @code{v} has not already been registered
         */
        void newVar(Var v, bool sign, bool dvar);

        /**
         * @brief Declare whether a variable should be eligible for selection in the decision
         * heuristic.
         * 
         * @param v The variable to set
         * @param b true iff v is a decision variable
         */
        void setDecisionVar(Var v, bool b); 

        /**
         * @brief Return the next decision variable.
         * 
         * @return the decision literal
         */
        Lit pickBranchLit(void);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // VARIABLE SELECTION HEURISTIC STATE MODIFICATION

        /**
         * @brief Rebuild the priority queue from scratch
         * 
         */
        void rebuildPriorityQueue(void);

        /**
         * @brief Set the variable activity
         * 
         * @param v the variable whose activity should be set
         * @param activityValue the activity to give the variable
         */
        void setActivity(Var v, double activityValue);

    #if BRANCHING_HEURISTIC == VSIDS
        // VSIDS

        /**
         * @brief Decay all variables with the specified factor. 
         * 
         * @note Implemented by increasing the 'bump' value instead.
         */
        void decayActivityVSIDS(void);

        /**
         * @brief Increase a variable with the current 'bump' value.
         * 
         * @param v The variable to bump
         */
        void bumpActivityVSIDS (Var v);
    #endif

        /**
         * @brief Update the position of a variable in the priority queue with respect to a
         * decreased activity
         * 
         * @param v the variable whose activity was decreased
         */
        void decreasePriorityQueue(Var v);

        /**
         * @brief Update the position of a variable in the priority queue with respect to an
         * increased activity
         * 
         * @param v the variable whose activity was increased
         */
        void increasePriorityQueue(Var v);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // POLARITY SELECTION HEURISTIC STATE MODIFICATION

        /**
         * @brief Declare which polarity the decision heuristic should use for a variable. Requires
         * mode 'polarity_user'.
         * 
         * @param v The variable to set
         * @param b The default polarity for v
         */
        void setPolarity   (Var v, bool b);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // EVENT LISTENERS

        /**
         * @brief Update data structures for branching heuristics upon receiving a new input clause
         * 
         * @param ps The literal that was assigned
         */
        void handleEventInputClause(const vec<Lit>& ps);

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
         * @param assignedAtLastLevel true iff l was assigned at the last decision level
         */
        void handleEventLitUnassigned(Lit l, uint64_t conflicts, bool assignedAtLastLevel);

        /**
         * @brief Update data structures for branching heuristics when unit propagation exits.
         * 
         * @param conflicts the new total number of conflicts
         * @param consistentState true iff propagation did not result in conflict
         */
        void handleEventPropagated(uint64_t conflicts, bool consistentState);

        /**
         * @brief Update data structures for branching heuristics upon conflict.
         * 
         * @param conflicts the new total number of conflicts
         */
        void handleEventConflicted(uint64_t conflicts);

        /**
         * @brief Update data structures for branching heuristics when a literal appears in the
         * conflict graph.
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
         * @param out_learnt the simplified learnt clause
         * @param seen true iff a variable was in the original learnt clause
         * 
         * @post @code{seen} is unchanged
         */
        void handleEventLearnedClause(const vec<Lit>& out_learnt, vec<bool>& seen);

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // HELPER FUNCTIONS

        /**
         * @brief Insert a variable in the decision order priority queue.
         * 
         * @param x the variable to insert
         */
        void insertVarOrder(Var x);

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // EVENT LISTENER HELPER FUNCTIONS

    #if BRANCHING_HEURISTIC == CHB
        /**
         * @brief Update data structures for CHB upon variable assignment.
         * 
         * @param l The literal that was assigned
         * @param conflicts The current number of conflicts in the solver
         */
        void handleEventLitAssignedCHB(Lit l, uint64_t conflicts);

        /**
         * @brief Update data structures for CHB upon variable unassignment.
         * 
         * @param l The literal that was unassigned
         * @param conflicts The current number of conflicts in the solver
         * @param assignedAtLastLevel true iff l was assigned at the last decision level
         */
        void handleEventLitUnassignedCHB(Lit l, uint64_t conflicts, bool assignedAtLastLevel);

    #elif BRANCHING_HEURISTIC == LRB
        /**
         * @brief Update data structures for LRB upon picking a branch literal using LRB
         * 
         * @param conflicts The current number of conflicts in the solver
         */
        void handleEventPickBranchLitLRB(uint64_t conflicts);

        /**
         * @brief Update data structures for LRB upon variable assignment.
         * 
         * @param l The literal that was assigned
         * @param conflicts The current number of conflicts in the solver
         */
        void handleEventLitAssignedLRB(Lit l, uint64_t conflicts);

        /**
         * @brief Update data structures for LRB upon variable unassignment.
         * 
         * @param l The literal that was unassigned
         * @param conflicts The current number of conflicts in the solver
         * @param assignedAtLastLevel true iff l was assigned at the last decision level
         */
        void handleEventLitUnassignedLRB(Lit l, uint64_t conflicts, bool assignedAtLastLevel);
    #endif
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE FUNCTIONS

    //////////////
    // ACCESSORS

    inline const vec<double>& BranchingHeuristicManager::getActivityVSIDS() const {
        return activity;
    }

    ///////////////////////
    // STATE MODIFICATION

    inline void BranchingHeuristicManager::newVar(Var v, bool sign, bool dvar) {
        // Decision variables
        decision.push();
        setDecisionVar(v, dvar);
        
        // VSIDS
        activity.push(rnd_init_act ? randomNumberGenerator.drand() * 0.00001 : 0);

    #if BRANCHING_HEURISTIC == CHB
        // CHB
        conflicted.push(0);
        picked.push(0);
        last_conflict.push(0);

    #elif BRANCHING_HEURISTIC == LRB
        // LRB
        conflicted.push(0);
        picked.push(0);
    #if ALMOST_CONFLICT
        almost_conflicted.push(0);
    #endif
    #if ANTI_EXPLORATION
        canceled.push(0);
    #endif
    #endif

        // Phase saving
        polarity.push(sign);

    #if PRIORITIZE_ER
        // Extension level branching
        degree.push(0);
        extensionLevel.push(0);
    #endif

    #if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
        // Statistics
        total_actual_rewards.push(0);
        total_actual_count.push(0);
    #endif
    }

    inline void BranchingHeuristicManager::setDecisionVar(Var v, bool b) { 
        if      ( b && !decision[v]) dec_vars++;
        else if (!b &&  decision[v]) dec_vars--;

        decision[v] = b;
        insertVarOrder(v);
    }

    ////////////////////////////////////////////////////
    // VARIABLE SELECTION HEURISTIC STATE MODIFICATION

    inline void BranchingHeuristicManager::setActivity(Var v, double activityValue) {
        const double oldActivity = activity[v];
        activity[v] = activityValue;
        if (activity[v] > oldActivity)
            increasePriorityQueue(v);
        else
            decreasePriorityQueue(v);
    }

#if BRANCHING_HEURISTIC == VSIDS
    inline void BranchingHeuristicManager::decayActivityVSIDS() {
        var_inc *= (1 / var_decay);
    }

    inline void BranchingHeuristicManager::bumpActivityVSIDS(Var v) {
        const double RESCALE_THRESHOLD = 1e100;

        if ((activity[v] += var_inc) > RESCALE_THRESHOLD) {
        #if PRIORITIZE_ER && defined(EXTLVL_ACTIVITY)
            // Clear extension level activity
            for (int i = 0; i < extensionLevelActivity.size(); i++) {
                extensionLevelActivity[i] = 0;
            }
        #endif
            // Rescale:
            for (int i = 0; i < assignmentTrail.nVars(); i++) {
                activity[i] /= RESCALE_THRESHOLD;
        #if PRIORITIZE_ER && defined(EXTLVL_ACTIVITY)
                extensionLevelActivity[extensionLevel[i]] += activity[i];
        #endif
            }
            var_inc /= RESCALE_THRESHOLD;
        }

        // Update variable in priority queue with respect to new activity
        increasePriorityQueue(v);
    }
#endif

    inline void BranchingHeuristicManager::decreasePriorityQueue(Var v) {
        // Update order_heap with respect to new activity
    #if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
        if (order_heap_extlvl.inHeap(v))
            order_heap_extlvl.increase(v);
        if (order_heap_degree.inHeap(v))
            order_heap_degree.increase(v);
    #else
        if (order_heap.inHeap(v))
            order_heap.increase(v);
    #endif
    }

    inline void BranchingHeuristicManager::increasePriorityQueue(Var v) {
        // Update order_heap with respect to new activity
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

    ////////////////////////////////////////////////////
    // POLARITY SELECTION HEURISTIC STATE MODIFICATION

    inline void BranchingHeuristicManager::setPolarity(Var v, bool b) {
        polarity[v] = b;
    }

    ////////////////////
    // EVENT LISTENERS

    inline void BranchingHeuristicManager::handleEventInputClause(const vec<Lit>& ps) {
    #if PRIORITIZE_ER
        for (int k = 0; k < ps.size(); k++)
            degree[var(ps[k])]++;
    #endif
    }

    inline void BranchingHeuristicManager::handleEventLitAssigned(Lit l, uint64_t conflicts) {
    #if BRANCHING_HEURISTIC == CHB
        handleEventLitAssignedCHB(l, conflicts);
    #elif BRANCHING_HEURISTIC == LRB
        handleEventLitAssignedLRB(l, conflicts);
    #endif
    }

    inline void BranchingHeuristicManager::handleEventLitUnassigned(
        Lit l,
        uint64_t conflicts,
        bool assignedAtLastLevel
    ) {
        const Var v = var(l);
    #if BRANCHING_HEURISTIC == CHB
        handleEventLitUnassignedCHB(l, conflicts, assignedAtLastLevel);
    #elif BRANCHING_HEURISTIC == LRB
        handleEventLitUnassignedLRB(l, conflicts, assignedAtLastLevel);
    #endif

        // Phase saving
        if (phase_saving == PhaseSavingLevel::FULL ||
            (phase_saving == PhaseSavingLevel::LIMITED) && assignedAtLastLevel
        ) {
            setPolarity(v, sign(l));
        }

        // Update priority queue
        insertVarOrder(v);
    }

    inline void BranchingHeuristicManager::handleEventPropagated(
        uint64_t conflicts,
        bool consistentState
    ) {
    #if BRANCHING_HEURISTIC == CHB
        const double multiplier = consistentState ? reward_multiplier : 1.0;
        for (int a = prev_trail_index; a < assignmentTrail.nAssigns(); a++) {
            Var v = var(assignmentTrail[a]);
            const uint64_t age = conflicts - last_conflict[v] + 1;
            const double reward = multiplier / age ;
            const double old_activity = activity[v];
            activity[v] = step_size * reward + ((1 - step_size) * old_activity);
            if (activity[v] > old_activity) increasePriorityQueue(v);
            else                            decreasePriorityQueue(v);
        }

        prev_trail_index = assignmentTrail.nAssigns();
    #endif
    }

    inline void BranchingHeuristicManager::handleEventConflicted(uint64_t conflicts) {
    #if BRANCHING_HEURISTIC == VSIDS
        decayActivityVSIDS();
    #elif BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
        if (step_size > min_step_size)
            step_size -= step_size_dec;
    #endif

    #ifdef POLARITY_VOTING
        // Count votes for the polarity that led to the conflict
        polarity_count.clear();
        for (int k = 0; k < group_polarity.size(); k++) polarity_count.push(0);
    #endif
    }

    inline void BranchingHeuristicManager::handleEventLitInConflictGraph(
        Lit q,
        uint64_t conflicts
    ) {
        const Var v = var(q);

        // Update data structures for variable selection heuristic
    #if BRANCHING_HEURISTIC == VSIDS
        bumpActivityVSIDS(v);
    #elif BRANCHING_HEURISTIC == CHB
        last_conflict[v] = conflicts;
        conflicted[v]++;
    #elif BRANCHING_HEURISTIC == LRB
        conflicted[v]++;
    #endif

        // Update data structures for polarity selection heuristic
    #ifdef POLARITY_VOTING
        // Count votes for the polarity that led to the conflict
        count[extensionLevel[v]] += sign(q) ? (+1) : (-1);
    #endif
    }

    inline void BranchingHeuristicManager::handleEventRestarted(uint64_t propagations) {

    }

    /////////////////////
    // HELPER FUNCTIONS

    inline void BranchingHeuristicManager::insertVarOrder(Var x) {
    #if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
        Heap<VarOrderLt>& order_heap = extensionLevel[x] ? order_heap_extlvl : order_heap_degree;
    #endif
        if (!order_heap.inHeap(x) && decision[x]) order_heap.insert(x);
    }

    ////////////////////////////////////
    // EVENT LISTENER HELPER FUNCTIONS

#if BRANCHING_HEURISTIC == CHB
    inline void BranchingHeuristicManager::handleEventLitAssignedCHB(Lit l, uint64_t conflicts) {
        const Var v = var(l);
        picked[v] = conflicts;
        conflicted[v] = 0;
    }

    inline void BranchingHeuristicManager::handleEventLitUnassignedCHB(
        Lit l,
        uint64_t conflicts,
        bool assignedAtLastLevel
    ) {
        const Var v = var(l);
        uint64_t age = conflicts - picked[v];
        if (age > 0) {
            const double reward = ((double) conflicted[v]) / ((double) age);
            const double old_activity = activity[v];
            activity[v] = step_size * reward + ((1 - step_size) * old_activity);
            if (activity[v] > old_activity) increasePriorityQueue(v);
            else                            decreasePriorityQueue(v);

            // Update statistics
            total_actual_rewards[v] += reward;
            total_actual_count[v] ++;
        }
    }

#elif BRANCHING_HEURISTIC == LRB
    inline void BranchingHeuristicManager::handleEventPickBranchLitLRB(uint64_t conflicts) {
    #if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
        auto& order_heap = order_heap_extlvl.empty() ? order_heap_degree : order_heap_extlvl;
    #endif
    #if ANTI_EXPLORATION
        Var next = order_heap[0];
        uint64_t age = conflicts - canceled[next];
        while (age > 0 && assignmentTrail.value(next) == l_Undef) {
            double decay = pow(0.95, age);
            activity[next] *= decay;
            if (order_heap.inHeap(next)) {
                order_heap.increase(next);
            }
            canceled[next] = conflicts;
            next = order_heap[0];
            age = conflicts - canceled[next];
        }
    #endif
    }

    inline void BranchingHeuristicManager::handleEventLitAssignedLRB(Lit l, uint64_t conflicts) {
        const Var v = var(l);
        picked[v] = conflicts;
    #if ANTI_EXPLORATION
        uint64_t age = conflicts - canceled[v];
        if (age > 0) {
            double decay = pow(0.95, age);
            activity[v] *= decay;
            decreasePriorityQueue(v);
        }
    #endif
        conflicted[v] = 0;
    #if ALMOST_CONFLICT
        almost_conflicted[v] = 0;
    #endif
    }

    inline void BranchingHeuristicManager::handleEventLitUnassignedLRB(
        Lit l,
        uint64_t conflicts,
        bool assignedAtLastLevel
    ) {
        const Var v = var(l);
        uint64_t age = conflicts - picked[v];
        if (age > 0) {
            const double reward = ((double) conflicted[v]) / ((double) age);
        #if ALMOST_CONFLICT
            const double adjusted_reward = ((double) (conflicted[v] + almost_conflicted[v])) / age;
        #else
            const double adjusted_reward = reward;
        #endif
            const double old_activity = activity[v];
            activity[v] = step_size * adjusted_reward + ((1 - step_size) * old_activity);
            if (activity[v] > old_activity) increasePriorityQueue(v);
            else                            decreasePriorityQueue(v);

            // Update statistics
            total_actual_rewards[v] += reward;
            total_actual_count[v] ++;
        }
    #if ANTI_EXPLORATION
        canceled[v] = conflicts;
    #endif
    }
#endif

} // namespace Minisat

#endif