/***********************************************************************[RestartHeuristicManager.h]
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

#ifndef Minisat_RestartHeuristicManager_h
#define Minisat_RestartHeuristicManager_h

#include <math.h>
#include "core/SolverTypes.h"

namespace Minisat {
    // Forward declarations
    class Solver;

    class RestartHeuristicManager {
    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // HELPER TYPES

        template<typename T>
        class MyQueue {
            int max_sz, q_sz;
            int ptr;
            int64_t sum;
            vec<T> q;
        public:
            MyQueue(int sz) : max_sz(sz), q_sz(0), ptr(0), sum(0) { assert(sz > 0); q.growTo(sz); }
            inline bool   full () const { return q_sz == max_sz; }
    #ifdef INT_QUEUE_AVG
            inline T      avg  () const { assert(full()); return sum / max_sz; }
    #else
            inline double avg  () const { assert(full()); return sum / (double) max_sz; }
    #endif
            inline void   clear()       { sum = 0; q_sz = 0; ptr = 0; }
            void push(T e) {
                if (q_sz < max_sz) q_sz++;
                else sum -= q[ptr];
                sum += e;
                q[ptr++] = e;
                if (ptr == max_sz) ptr = 0;
            }
        };
    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // MEMBER VARIABLES

        const bool& VSIDS;

        bool cached;

        float global_lbd_sum;
        
        /// @brief For computing moving averages of recent LBD values.
        MyQueue<int> lbd_queue;

        int conflictBudget;


    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////////
        // PARAMETERS

        /// @brief The initial restart limit. (default 100)
        int restart_first;

        /// @brief The factor by which the restart limit is multiplied in each restart. (default 1.5)
        double restart_inc;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATISTICS

        /// @brief Current number of restarts
        mutable int curr_restarts;

        /// @brief The total number of conflicts encountered during search while using VSIDS.
        uint64_t conflicts_VSIDS;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CONSTRUCTORS

        RestartHeuristicManager(Solver& s);
        ~RestartHeuristicManager() = default;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // ACCESSORS

        bool shouldRestart(void);

        int getConflictBudget(void) const;

        int getRestartConflicts(void) const;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // EVENT HANDLERS

        void handleEventRestarted(int nof_conflicts);

        void handleEventLearntClause(int lbd);
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE METHODS

    //////////////
    // ACCESSORS

    inline bool RestartHeuristicManager::shouldRestart() {
        if (!VSIDS)
            return conflictBudget <= 0;
        else if (!cached){
            cached = true;
            return lbd_queue.full() && (lbd_queue.avg() * 0.8 > global_lbd_sum / conflicts_VSIDS);
        }

        return false;
    }

    inline int RestartHeuristicManager::getConflictBudget(void) const {
        return conflictBudget;
    }

    /*
    Finite subsequences of the Luby-sequence:

    0: 1
    1: 1 1 2
    2: 1 1 2 1 1 2 4
    3: 1 1 2 1 1 2 4 1 1 2 1 1 2 4 8
    ...


    */

    inline static double luby(double y, int x) {
        // Find the finite subsequence that contains index 'x', and the
        // size of that subsequence:
        int size, seq;
        for (size = 1, seq = 0; size < x + 1; seq++, size = 2 * size + 1);

        while (size - 1 != x) {
            size = (size - 1) >> 1;
            seq--;
            x = x % size;
        }

        return pow(y, seq);
    }

    inline int RestartHeuristicManager::getRestartConflicts(void) const {
        if (VSIDS) {
            return INT32_MAX;
        } else {
            return luby(restart_inc, curr_restarts++) * restart_first;
        }
    }

    ///////////////////
    // EVENT HANDLERS

    inline void RestartHeuristicManager::handleEventRestarted(int nof_conflicts) {
        conflictBudget = nof_conflicts;
        lbd_queue.clear();
        cached = false;
    }

    inline void RestartHeuristicManager::handleEventLearntClause(int lbd) {
        conflictBudget--;
        if (VSIDS) {
            cached = false;
            conflicts_VSIDS++;
            lbd_queue.push(lbd);
            global_lbd_sum += (lbd > 50 ? 50 : lbd);
        }
    }
}

#endif