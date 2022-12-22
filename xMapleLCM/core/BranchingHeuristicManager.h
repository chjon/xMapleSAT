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

#define ANTI_EXPLORATION

#include <math.h>
#include "core/SolverTypes.h"
#include "core/SolverERTypes.h"
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
        vec<double>
            activity_VSIDS,
            activity_CHB,
            activity_distance;

        // A priority queue of variables ordered with respect to the variable activity.
        Heap<VarOrderLt>
            order_heap_VSIDS,
            order_heap_CHB,
            order_heap_distance;

        vec<char> decision; // Declares whether a variable is eligible for selection in the decision heuristic.
        
        //////////////////////////
        // Heuristic configuration

        // VSIDS
        int    timer;   // Frequency for updating var_decay, measured in conflicts 
        double var_inc; // Amount to bump next variable with.
        double var_decay;
        vec<Lit> conflictLits; // Literals that participate in the conflict graph

        // CHB
        double step_size;
        double step_size_dec;
        double min_step_size;
        vec<uint32_t> conflicted;
        vec<uint32_t> almost_conflicted;
        vec<uint32_t> picked;
#ifdef ANTI_EXPLORATION
        vec<uint32_t> canceled;
#endif

        // Distance
        vec<Lit>    involved_lits;
        vec<int>    pathCs;
        vec<double> var_iLevel;
        vec<double> var_iLevel_tmp;
        double      var_iLevel_inc;
        double      my_var_decay;

        // Random
        double random_var_freq;
        bool   rnd_pol;      // Use random polarities for branching heuristics.
        bool   rnd_init_act; // Initialize variable activities with a small random value.

        // Phase saving
        int phase_saving;   // Controls the level of phase saving (0=none, 1=limited, 2=full).
        vec<char> polarity; // The preferred polarity of each variable.

        // Current branching heuristic
        uint64_t VSIDS_props_limit;
        uint64_t prev_props;
        bool switch_mode;
        bool m_VSIDS;
        bool m_DISTANCE;

    public:
        ////////////////
        // STATISTICS //
        ////////////////

        uint64_t dec_vars     ; // Total number of decision variables
        uint64_t decisions    ; // Total number of decisions
        uint64_t rnd_decisions; // Total number of random decisions

        uint64_t conflicts_VSIDS;

    protected:
        ///////////////////////
        // SOLVER REFERENCES //
        ///////////////////////

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

        // VSIDS

        void decayActivityVSIDS();                   // Decay all variables with the specified factor. Implemented by increasing the 'bump' value instead.
        void bumpActivityVSIDS (Var v, double mult); // Increase a variable with the current 'bump' value.

        // CHB

        void clearCHBStats(); // Clear CHB statistics

        // Distance
        void bumpActivityDistance(vec<Lit>& involved_lits, int max_level);
        void collectFirstUIP(CRef confl);

        // Extended resolution

        /**
         * @brief Prioritize branching on a given set of variables
         * 
         * @param defs Extension variable definitions -- the extension variables will be prioritized
         * @param scaleFactor The fraction of the current top activity to assign to the new variables
         */
        void prioritize(const std::vector<ExtDef>& defs, double scaleFactor);

        // Choose branching heuristic
        void setBranchingHeuristic(BranchingHeuristic branchingHeuristic);

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
         */
        void handleEventLitInConflictGraph(Lit l);

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
        BranchingHeuristic currentBranchingHeuristic() const;
        bool distanceBranchingEnabled() const;

    private:
        /////////////////////////////////////
        // HELPER FUNCTIONS FOR HEURISTICS //
        /////////////////////////////////////
        void collectFirstUIPConflictClause(int& minLevel, CRef confl);
        void collectFirstUIPReasonClause  (int& minLevel, CRef confl, int reasonVarLevel);
    };

    inline void BranchingHeuristicManager::insertVarOrder(Var x) {
        //    Heap<VarOrderLt>& order_heap = m_VSIDS ? order_heap_VSIDS : order_heap_CHB;
        auto& order_heap = m_DISTANCE ? order_heap_distance : ((!m_VSIDS) ? order_heap_CHB:order_heap_VSIDS);
        if (!order_heap.inHeap(x) && decision[x]) order_heap.insert(x);
    }

    inline void BranchingHeuristicManager::decayActivityVSIDS() { if (m_VSIDS) var_inc *= (1 / var_decay); }

    inline void BranchingHeuristicManager::bumpActivityVSIDS(Var v, double mult) {
        if ((activity_VSIDS[v] += var_inc * mult) > 1e100) {
            // Rescale:
            for (int i = 0; i < variableDatabase.nVars(); i++)
                activity_VSIDS[i] *= 1e-100;
            var_inc *= 1e-100;
        }

        // Update order_heap with respect to new activity:
        if (order_heap_VSIDS.inHeap(v)) order_heap_VSIDS.decrease(v);
        unitPropagator.increasePriority(v);
    }

    inline void BranchingHeuristicManager::clearCHBStats() {
        // Instead of clearing, set vectors to 0
        for (int i = 0; i < variableDatabase.nVars(); i++) {
            picked           [i] = 0;
            conflicted       [i] = 0;
            almost_conflicted[i] = 0;
#ifdef ANTI_EXPLORATION
            canceled         [i] = 0;
#endif
        }
    }

    inline void BranchingHeuristicManager::bumpActivityDistance(vec<Lit>& involved_lits, int max_level) {
        double inc=var_iLevel_inc;
        vec<int> level_incs; level_incs.clear();
        for(int i = 0; i < max_level; i++){
            level_incs.push(inc);
            inc = inc / my_var_decay;
        }

        for (int i = 0; i < involved_lits.size(); i++) {
            Var v = var(involved_lits[i]);
            // double old_act=activity_distance[v];
            // activity_distance[v] +=var_iLevel_inc * var_iLevel_tmp[v];
            activity_distance[v] += var_iLevel_tmp[v]*level_incs[var_iLevel_tmp[v]-1];

            if (activity_distance[v] > 1e100) {
                for (int vv = 0; vv < variableDatabase.nVars(); vv++)
                    activity_distance[vv] *= 1e-100;
                
                var_iLevel_inc *= 1e-100;
                for (int j = 0; j < max_level; j++)
                    level_incs[j] *= 1e-100;
            }
            if (order_heap_distance.inHeap(v))
                order_heap_distance.decrease(v);

            // var_iLevel_inc *= (1 / my_var_decay);
        }
        var_iLevel_inc=level_incs[level_incs.size()-1];

    }

    inline void BranchingHeuristicManager::setBranchingHeuristic(BranchingHeuristic branchingHeuristic) {
        switch (branchingHeuristic) {
            case BranchingHeuristic::VSIDS: { m_VSIDS = true;  } break;
            case BranchingHeuristic::CHB:   { m_VSIDS = false; } break;
            default: { assert(0); break; } // Other heuristics are not supported
        }
    }

    inline void BranchingHeuristicManager::setPolarity   (Var v, bool b) { polarity[v] = b; }
    inline void BranchingHeuristicManager::setDecisionVar(Var v, bool b) {
        if      ( b && !decision[v]) dec_vars++;
        else if (!b &&  decision[v]) dec_vars--;

        decision[v] = b;
        if (b && !order_heap_CHB     .inHeap(v)) order_heap_CHB     .insert(v);
        if (b && !order_heap_VSIDS   .inHeap(v)) order_heap_VSIDS   .insert(v);
        if (b && !order_heap_distance.inHeap(v)) order_heap_distance.insert(v);
    }

    inline void BranchingHeuristicManager::handleEventLitAssigned(Lit l, uint64_t conflicts) {
        const Var v = var(l);

        // CHB
        if (!m_VSIDS){
            picked[v] = conflicts;
            conflicted[v] = 0;
            almost_conflicted[v] = 0;
    #ifdef ANTI_EXPLORATION
            uint32_t age = conflicts - canceled[v];
            if (age > 0) {
                double decay = pow(0.95, age);
                activity_CHB[v] *= decay;
                if (order_heap_CHB.inHeap(v))
                    order_heap_CHB.increase(v);
            }
    #endif
        }
    }

    inline void BranchingHeuristicManager::handleEventLitUnassigned(Lit l, uint64_t conflicts, bool assignedAtLastLevel) {
        const Var v = var(l);

        // CHB
        if (!m_VSIDS){
            uint32_t age = conflicts - picked[v];
            if (age > 0) {
                double adjusted_reward = ((double) (conflicted[v] + almost_conflicted[v])) / ((double) age);
                double old_activity = activity_CHB[v];
                activity_CHB[v] = step_size * adjusted_reward + ((1 - step_size) * old_activity);
                if (order_heap_CHB.inHeap(v)) {
                    if (activity_CHB[v] > old_activity) order_heap_CHB.decrease(v);
                    else                                order_heap_CHB.increase(v);
                }
            }
#ifdef ANTI_EXPLORATION
            canceled[v] = conflicts;
#endif
        }

        // Phase saving
        if (phase_saving > 1 || (phase_saving == 1) && assignedAtLastLevel)
            setPolarity(v, sign(l));

        // Update priority queue
        insertVarOrder(v);
    }

    inline void BranchingHeuristicManager::handleEventConflicted(uint64_t conflicts) {
        if (m_VSIDS) {
            if (--timer == 0 && var_decay < 0.95) {
                timer = 5000;
                var_decay += 0.01;
            }

            conflicts_VSIDS++;
        } else{
            if (step_size > min_step_size)
                step_size -= step_size_dec;
        }

        // Enable distance branching
        if (conflicts > 50000) m_DISTANCE = 0;
        else m_DISTANCE = 1;
    }

    inline void BranchingHeuristicManager::handleEventLitInConflictGraph(Lit q) {
        if (m_VSIDS) {
            bumpActivityVSIDS(var(q), .5);
            conflictLits.push(q);
        } else {
            conflicted[var(q)]++;
        }
    }

    inline void BranchingHeuristicManager::handleEventRestarted(uint64_t propagations) {
        // Switch branching heuristics
        if (switch_mode) { 
            switch_mode = false;
            m_VSIDS = !m_VSIDS;
            if (m_VSIDS) printf("c Switched to VSIDS.\n");
            else         printf("c Switched to LRB.\n");
            fflush(stdout);
            clearCHBStats();
        }

        // Schedule next switch
        if (propagations - prev_props > VSIDS_props_limit){
            prev_props = propagations;
            VSIDS_props_limit = VSIDS_props_limit + VSIDS_props_limit / 10;
            switch_mode = true;
        }
    }

    inline const vec<double>& BranchingHeuristicManager::getActivityVSIDS() const { return activity_VSIDS; }
    inline BranchingHeuristic BranchingHeuristicManager::currentBranchingHeuristic() const { return m_VSIDS ? BranchingHeuristic::VSIDS : BranchingHeuristic::CHB; }
    inline bool BranchingHeuristicManager::distanceBranchingEnabled() const { return m_DISTANCE; }
} // namespace Minisat

#endif