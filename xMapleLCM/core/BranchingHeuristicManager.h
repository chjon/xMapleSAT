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

#define ANTI_EXPLORATION

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

        /// @brief The CHB activity of a variable.
        vec<double> activity_CHB;

        /// @brief The VSIDS activity of a variable.
        vec<double> activity_VSIDS;

        /// @brief The distance activity of a variable.
        vec<double> activity_distance;

        /// @brief A priority queue of variables ordered with respect to the variable activity.
        Heap< VarOrderLt<double> > order_heap_CHB;
        Heap< VarOrderLt<double> > order_heap_VSIDS;
        Heap< VarOrderLt<double> > order_heap_distance;

        //////////
        // VSIDS


        /// @brief Amount by which to bump variable activity
        double var_inc;

        /// @brief Amount by which to decay variable activities
        double var_decay;

        int timer;

        ////////
        // CHB

        /// @brief Step size for computing CHB activity
        double step_size;

        double step_size_dec;
        double min_step_size;

        ////////
        // LRB

    public:
        /// @brief The currently selected heuristic - true for VSIDS, false for CHB
        bool VSIDS;

    protected:
        /// @brief The number of times a variable appears in a conflict graph
        vec<uint32_t> conflicted;

        /// @brief The total number of conflicts seen by the solver before the variable was
        /// assigned
        vec<uint32_t> picked;

        /// @brief The number of times a variable appeared directly before the cut in the conflict
        /// graph for a learnt clause
        vec<uint32_t> almost_conflicted;

    #ifdef ANTI_EXPLORATION
        /// @brief The total number of conflicts seen by the solver before the variable was
        /// selected as a decision literal or unassigned
        vec<uint32_t> canceled;
    #endif

        /////////////
        // DISTANCE

    public:
        bool DISTANCE;

    protected:
        double var_iLevel_inc;

        double my_var_decay;

        ///////////////////////////////////////////////////////////////////////////////////////////
        // POLARITY SELECTION HEURISTIC MEMBER VARIABLES

        /// @brief The preferred polarity of each variable.
        vec<char> polarity;

        ///////////////////
        // MODE SWITCHING

        bool switchMode;

        uint64_t prevPropagations;
        
    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // PARAMETERS

        ///////////
        // Random

        /// @brief Initialize variable activities with a small random value
        bool rnd_init_act;

        /////////////////
        // Phase saving

        /// @brief Controls the level of phase saving (0=none, 1=limited, 2=full).
        PhaseSavingLevel phase_saving;

        /// @brief Specifies the number of propagations after which the solver switches between LRB and VSIDS
        uint64_t VSIDS_props_limit;

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
         * @brief Get the current activity array
         * 
         * @return the current activity array 
         */
        const vec<double>& getActivity(void) const;

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

        /**
         * @brief Check whether to switch the current branching heuristic
         * 
         * @param propagations 
         */
        void checkSwitchHeuristic(uint64_t propagations);

        /**
         * @brief Switch the current branching heuristic from VSIDS to LRB and vice versa.
         * 
         */
        void switchHeuristic(void);

        void updateActivityDistance(
            const vec<Var>& involvedVars,
            const vec<int>& varDist,
            int max_distance
        );

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
         * @param confl the conflicting clause
         * @param conflicts the new total number of conflicts
         */
        void handleEventConflicted(CRef confl, uint64_t conflicts);

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
         * @param backtrackLevel the backtrack level for the learnt clause
         * 
         * @post @code{seen} is unchanged
         */
        void handleEventLearnedClause(const vec<Lit>& out_learnt, vec<bool>& seen, int backtrackLevel);

        void handleEventPickBranchLit(uint64_t conflicts);

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

        void varDecayActivity(void);

        void varBumpActivity(Var v, double mult);
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE FUNCTIONS

    //////////////
    // ACCESSORS

    inline const vec<double>& BranchingHeuristicManager::getActivity() const {
        if (DISTANCE) return activity_distance;
        if (VSIDS)    return activity_VSIDS;
        return activity_CHB;
    }

    ///////////////////////
    // STATE MODIFICATION

    inline void BranchingHeuristicManager::newVar(Var v, bool sign, bool dvar) {
        // Decision variables
        decision.push();
        setDecisionVar(v, dvar);
        
        // Branching heuristics
        activity_CHB     .push(0);
        activity_VSIDS   .push(rnd_init_act ? randomNumberGenerator.drand() * 0.00001 : 0);
        activity_distance.push(0);

        picked.push(0);
        conflicted.push(0);
        almost_conflicted.push(0);
    #ifdef ANTI_EXPLORATION
        canceled.push(0);
    #endif

        // Phase saving
        polarity.push(sign);
    }

    inline void BranchingHeuristicManager::setDecisionVar(Var v, bool b) { 
        if      ( b && !decision[v]) dec_vars++;
        else if (!b &&  decision[v]) dec_vars--;

        decision[v] = b;
        if (b && !order_heap_CHB.inHeap(v)) {
            order_heap_CHB     .insert(v);
            order_heap_VSIDS   .insert(v);
            order_heap_distance.insert(v);
        }
    }

    ////////////////////////////////////////////////////
    // VARIABLE SELECTION HEURISTIC STATE MODIFICATION

    inline void BranchingHeuristicManager::setActivity(Var v, double activityValue) {

    }

    inline void BranchingHeuristicManager::decreasePriorityQueue(Var v) {

    }

    inline void BranchingHeuristicManager::increasePriorityQueue(Var v) {

    }

    ////////////////////////////////////////////////////
    // POLARITY SELECTION HEURISTIC STATE MODIFICATION

    inline void BranchingHeuristicManager::setPolarity(Var v, bool b) {
        polarity[v] = b;
    }

    ////////////////////
    // EVENT LISTENERS

    inline void BranchingHeuristicManager::handleEventInputClause(const vec<Lit>& ps) {

    }

    inline void BranchingHeuristicManager::handleEventLitAssigned(Lit l, uint64_t conflicts) {
        // Nothing to do for VSIDS
        if (VSIDS) return;

        // Update CHB data structures
        Var x = var(l);
        picked[x] = conflicts;
        conflicted[x] = 0;
        almost_conflicted[x] = 0;
    #ifdef ANTI_EXPLORATION
        uint32_t age = conflicts - canceled[x];
        if (age > 0){
            double decay = pow(0.95, age);
            activity_CHB[x] *= decay;
            if (order_heap_CHB.inHeap(x))
                order_heap_CHB.increase(x);
        }
    #endif
    }

    inline void BranchingHeuristicManager::handleEventLitUnassigned(
        Lit l,
        uint64_t conflicts,
        bool assignedAtLastLevel
    ) {
        // Update CHB data structures
        Var x = var(l);
        if (!VSIDS) {
            uint32_t age = conflicts - picked[x];
            if (age > 0){
                double adjusted_reward = ((double) (conflicted[x] + almost_conflicted[x])) / ((double) age);
                double old_activity = activity_CHB[x];
                activity_CHB[x] = step_size * adjusted_reward + ((1 - step_size) * old_activity);
                if (order_heap_CHB.inHeap(x)){
                    if (activity_CHB[x] > old_activity)
                        order_heap_CHB.decrease(x);
                    else
                        order_heap_CHB.increase(x);
                }
            }
#ifdef ANTI_EXPLORATION
            canceled[x] = conflicts;
#endif
        }
        
        // Phase saving
        if (phase_saving == PhaseSavingLevel::FULL ||
            (phase_saving == PhaseSavingLevel::LIMITED) && assignedAtLastLevel
        ) {
            setPolarity(x, sign(l));
        }
        
        // Update priority queue
        insertVarOrder(x);
    }

    inline void BranchingHeuristicManager::handleEventPropagated(
        uint64_t conflicts,
        bool consistentState
    ) {

    }

    inline void BranchingHeuristicManager::handleEventLitInConflictGraph(
        Lit q,
        uint64_t conflicts
    ) {
        const Var v = var(q);
        if (VSIDS){
            varBumpActivity(v, .5);
            toClear.push(v);
        } else {
            conflicted[v]++;
        }
    }

    inline void BranchingHeuristicManager::handleEventRestarted(uint64_t propagations) {

    }

    inline void BranchingHeuristicManager::handleEventPickBranchLit(uint64_t conflicts) {
    #ifdef ANTI_EXPLORATION
        if (!VSIDS){
            Var v = order_heap_CHB[0];
            uint32_t age = conflicts - canceled[v];
            while (age > 0){
                double decay = pow(0.95, age);
                activity_CHB[v] *= decay;
                if (order_heap_CHB.inHeap(v))
                    order_heap_CHB.increase(v);
                canceled[v] = conflicts;
                v = order_heap_CHB[0];
                age = conflicts - canceled[v];
            }
        }
    #endif
    }

    /////////////////////
    // HELPER FUNCTIONS

    inline void BranchingHeuristicManager::insertVarOrder(Var x) {
        auto& order_heap = DISTANCE
            ? order_heap_distance
            : (VSIDS ? order_heap_VSIDS : order_heap_CHB);
        
        if (!order_heap.inHeap(x) && decision[x])
            order_heap.insert(x); 
    }

    inline void BranchingHeuristicManager::rebuildPriorityQueue(void) {
        vec<Var> vs;
        for (Var v = 0; v < assignmentTrail.nVars(); v++)
            if (decision[v] && assignmentTrail.value(v) == l_Undef)
                vs.push(v);

        order_heap_CHB     .build(vs);
        order_heap_VSIDS   .build(vs);
        order_heap_distance.build(vs);
    }

    ////////////////////////////////////
    // EVENT LISTENER HELPER FUNCTIONS

    inline void BranchingHeuristicManager::varDecayActivity() {
        var_inc *= (1 / var_decay);
    }

    inline void BranchingHeuristicManager::varBumpActivity(Var v, double mult) {
        constexpr double RESCALE_THRESHOLD = 1e100;
        if ((activity_VSIDS[v] += var_inc * mult) > RESCALE_THRESHOLD) {
            // Rescale:
            for (int i = 0; i < assignmentTrail.nVars(); i++)
                activity_VSIDS[i] /= RESCALE_THRESHOLD;
            var_inc /= RESCALE_THRESHOLD;
        }

        // Update order_heap with respect to new activity:
        if (order_heap_VSIDS.inHeap(v)) order_heap_VSIDS.decrease(v);
    }

} // namespace Minisat

#endif