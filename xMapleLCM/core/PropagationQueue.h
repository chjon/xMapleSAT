/******************************************************************************[PropagationQueue.h]
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

#ifndef Minisat_PropagationQueue_h
#define Minisat_PropagationQueue_h

// Define BCP-prioritization mode
#define BCP_PRIORITY_IMMEDIATE    0 // Immediate propagation
#define BCP_PRIORITY_DELAYED      1 // Delayed propagation
#define BCP_PRIORITY_OUT_OF_ORDER 3 // Out-of-order propagation

// Select initial prioritization mode
#ifndef BCP_PRIORITY_MODE
    #define BCP_PRIORITY_MODE BCP_PRIORITY_IMMEDIATE
#endif

// Define BCP-prioritization heuristic
#define BCP_PRIORITY_ACTIVITY   0 // Prioritize high-activity literals
#define BCP_PRIORITY_MAX_ON_MIN 1 // Prioritize variables that appear a lot on minimum-sized clauses

// Select prioritization heuristic
#ifndef BCP_PRIORITY_HEURISTIC
    #define BCP_PRIORITY_HEURISTIC BCP_PRIORITY_ACTIVITY
#endif

// Define BCP switching mode
#ifndef ENABLE_PRIORITY_BCP_RL
    #define ENABLE_PRIORITY_BCP_RL false
#endif

#include <limits.h>
#include "core/AssignmentTrail.h"
#include "core/BCPRLManager.h"
#include "core/SolverTypes.h"
#include "mtl/Heap.h"

namespace Minisat {
    // Forward declarations
    class Solver;

    /**
     * @brief This class represents the BCP propagation queue.
     * @note Implementation avoids polymorphism due to performance cost.
     */
    class PropagationQueue {
    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // HELPER TYPES

        /**
         * @brief Comparator for BCP priority queue
         */
        template<class T>
        struct LitOrderLt {
            const vec<T>*  activity;
            bool operator () (Var x, Var y) const {
                const vec<T>& act = *activity;
                x >>= 1; y >>= 1;
                return act[x] > act[y];
            }
            LitOrderLt(const vec<T>& act)
                : activity(&act)
            {}
        };

        struct OccurrenceCounter {
            int clauseSize;
            int count;

            bool operator < (const OccurrenceCounter other) const {
                if (clauseSize != other.clauseSize) return clauseSize < other.clauseSize;
                return count > other.count;
            }
            bool operator > (const OccurrenceCounter other) const {
                return !operator<(other);
            }
        };

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // SOLVER REFERENCES

        AssignmentTrail& assignmentTrail;
        BCPRLManager     bcprlManager;

    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // HEURISTICS

        // Maximum-Occurrences-on-Minimum-sized-Clauses heuristic
        vec<OccurrenceCounter> occurrences;

    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // MEMBER VARIABLES

        /// @brief The head of the propagation queue as an index into the assignment trail
        int qhead;

        /// @brief A view-only reference to the assignment trail
        const vec<Lit>& queue;

    #if BCP_PRIORITY_HEURISTIC == BCP_PRIORITY_MAX_ON_MIN
        /// @brief The priority queue for selecting variables to propagate during BCP
        Heap< LitOrderLt<OccurrenceCounter> > order_heap;
    #else // if BCP_PRIORITY_HEURISTIC == BCP_PRIORITY_ACTIVITY
        /// @brief The priority queue for selecting variables to propagate during BCP
        Heap< LitOrderLt<double> > order_heap;
    #endif

        /// @brief The polarities of queued variables
        vec<lbool> soft_assigns;

        /// @brief The reason clauses for queuing variables for BCP
        vec<CRef> reasons;

    public:
        /// @brief The current variant of BCP to use
        BCPMode current_bcpmode;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CONSTRUCTORS

        /**
         * @brief Construct a new PropagationQueue object
         * 
         * @param s Reference to main solver object
         */
        PropagationQueue(Solver& s);

        /**
         * @brief Destroy the PropagationQueue object
         * 
         */
        ~PropagationQueue() = default;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATE MODIFICATION

        /**
         * @brief Set up internal data structures for a new variable
         * 
         * @param v the variable to register
         */
        void newVar(Var v);

        /**
         * @brief Add a literal to the propagation queue
         * 
         * @param p the literal to add to the queue
         * @param from the reason for the literal
         * @return false if the negated version of the literal is already in the queue, true otherwise
         */
        bool enqueue(Lit p, CRef from = CRef_Undef);

        template <BCPMode bcpmode>
        bool enqueue(Lit p, CRef from = CRef_Undef);

        /**
         * @brief Add a literal to the propagation queue without notifying event listeners
         * 
         * @param p the literal to add to the queue
         * @param from the reason for the literal
         * @return false if the negated version of the literal is already in the queue, true otherwise
         */
        bool simpleEnqueue(Lit p, CRef from = CRef_Undef);

        template <BCPMode bcpmode>
        bool simpleEnqueue(Lit p, CRef from = CRef_Undef);

        /**
         * @brief Add a set of literals to the propagation queue
         * 
         * @param trail the current assignment trail
         * @param levelHead the start index from which to add literals to the queue
         * 
         * @note literals are added from @code{trail[levelHead]} up to @code{trail.last()}
         */
        void batchEnqueue(vec<Lit>& trail, int levelHead);

        /**
         * @brief Get the next literal for propagation
         * 
         * @tparam simple: true to skip updating stats and notifying event listeners
         * @return the literal to propagate
         */
        template <BCPMode bcpmode, bool simple>
        Lit getNext();

        /**
         * @brief Clear the propagation queue
         * 
         */
        template <BCPMode bcpmode>
        void clear();

        /**
         * @brief Set the activity metric to use for priority BCP
         * 
         * @param activity the activities of each variable
         */
        void prioritizeByActivity(const vec<double>& activity);
    
    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // EVENT HANDLERS

        void handleEventNewClause(const Clause& c);

        void handleEventRestarted(const BCPRLStats& stats);

        void handleEventLearntClause(uint64_t lbd);

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // HELPER FUNCTIONS

        /**
         * @brief Add a literal to the propagation queue
         * 
         * @tparam simple: true to skip updating stats and notifying event listeners
         * @param p the literal to add to the queue
         * @param from the reason for the literal
         * @return false if the negated version of the literal is already in the queue, true otherwise
         */
        template <BCPMode bcpmode, bool simple>
        bool genericEnqueue(Lit p, CRef from = CRef_Undef);
    };

    // Explicitly instantiate required templates
    template class PropagationQueue::LitOrderLt<double>;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE FUNCTIONS

    ///////////////////////
    // STATE MODIFICATION

    inline void PropagationQueue::newVar(Var v) {
        soft_assigns.push(l_Undef);
        reasons.push(CRef_Undef);
        occurrences.push(OccurrenceCounter{INT_MAX, 0});
    }

    template <BCPMode bcpmode>
    inline bool PropagationQueue::enqueue(Lit p, CRef from) {
        return genericEnqueue<bcpmode, false>(p, from);
    }

    inline bool PropagationQueue::enqueue(Lit p, CRef from) {
        switch (current_bcpmode) {
            default:
            case BCPMode::IMMEDIATE : return enqueue<BCPMode::IMMEDIATE >(p, from);
            case BCPMode::DELAYED   : return enqueue<BCPMode::DELAYED   >(p, from);
            case BCPMode::OUTOFORDER: return enqueue<BCPMode::OUTOFORDER>(p, from);
        }
    }

    template <BCPMode bcpmode>
    inline bool PropagationQueue::simpleEnqueue(Lit p, CRef from) {
        return genericEnqueue<bcpmode, true>(p, from);
    }

    inline bool PropagationQueue::simpleEnqueue(Lit p, CRef from) {
        switch (current_bcpmode) {
            default:
            case BCPMode::IMMEDIATE : return simpleEnqueue<BCPMode::IMMEDIATE >(p, from);
            case BCPMode::DELAYED   : return simpleEnqueue<BCPMode::DELAYED   >(p, from);
            case BCPMode::OUTOFORDER: return simpleEnqueue<BCPMode::OUTOFORDER>(p, from);
        }
    }

    inline void PropagationQueue::batchEnqueue(vec<Lit>& trail, int levelHead) {
        qhead = levelHead;
    }

    template <BCPMode bcpmode, bool simple>
    inline Lit PropagationQueue::getNext() {
        switch (bcpmode) {
            default:
            case BCPMode::IMMEDIATE: {
                return qhead == queue.size() ? lit_Undef : queue[qhead++];
            }
            case BCPMode::DELAYED: {
                // Propagate along the trail first
                if (qhead < queue.size()) return queue[qhead++];

                // Propagate according to priority queue afterward
                if (!order_heap.empty()) {
                    qhead++;
                    Lit p = Lit{order_heap.removeMin()};

                    // Note: Variable is assigned here because it is not assigned in @code{enqueue()}
                    if (simple) assignmentTrail.simpleAssign(p, reasons[var(p)]);
                    else        assignmentTrail.assign(p, reasons[var(p)]);
                    soft_assigns[var(p)] = l_Undef;
                    return p;
                }
                return lit_Undef;
            }
            case BCPMode::OUTOFORDER: {
                // Propagate along the trail first
                if (qhead < queue.size()) return queue[qhead++];

                // Check if propagation is done
                if (order_heap.empty()) return lit_Undef;

                // Propagate according to priority queue afterward
                qhead++;
                return Lit{order_heap.removeMin()};
            }
        }
    }

    template<>
    inline void PropagationQueue::clear<BCPMode::IMMEDIATE>() {
        qhead = queue.size();
    }

    template<>
    inline void PropagationQueue::clear<BCPMode::DELAYED>() {
        qhead = queue.size();
        for (int k = 0; k < order_heap.size(); k++)
            soft_assigns[order_heap[k] >> 1] = l_Undef;
        order_heap.clear();
    }

    template<>
    inline void PropagationQueue::clear<BCPMode::OUTOFORDER>() {
        order_heap.clear();
    }

    inline void PropagationQueue::prioritizeByActivity(const vec<double>& activity) {
    #if BCP_PRIORITY_MODE != BCP_PRIORITY_IMMEDIATE
    #if BCP_PRIORITY_HEURISTIC == BCP_PRIORITY_ACTIVITY
        order_heap.setComp(LitOrderLt<double>(activity));
    #endif
    #endif
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // HELPER FUNCTIONS

    template<BCPMode bcpmode, bool simple>
    inline bool PropagationQueue::genericEnqueue(Lit p, CRef from) {
        switch (bcpmode) {
            default:
            case BCPMode::IMMEDIATE: {
                if (simple) assignmentTrail.simpleAssign(p, from);
                else        assignmentTrail.assign(p, from);
            } break;
            case BCPMode::DELAYED: {
                if ((soft_assigns[var(p)] ^ sign(p)) == l_False) {
                    // Ensure conflicting literal is on the trail
                    if (assignmentTrail.value(p) == l_Undef) {
                        if (simple) assignmentTrail.simpleAssign(~p, reasons[var(p)]);
                        else        assignmentTrail.assign(~p, reasons[var(p)]);
                    }

                    return false;
                } else if ((soft_assigns[var(p)] ^ sign(p)) == l_Undef) {
                    // Note: actual variable assignment is delayed until the literal is popped by @code{getNext()}
                    order_heap.insert(p.x);
                    soft_assigns[var(p)] = lbool(!sign(p));
                    reasons[var(p)] = from;
                }
            } break;
            case BCPMode::OUTOFORDER: {
                order_heap.insert(p.x);
                if (simple) assignmentTrail.simpleAssign(p, from);
                else        assignmentTrail.assign(p, from);
            } break;
        }

        return true;
    }

    inline void PropagationQueue::handleEventNewClause(const Clause& c) {
    #if BCP_PRIORITY_MODE != BCP_PRIORITY_IMMEDIATE
        for (int i = 0; i < c.size(); i++) {
            Lit l = c[i];
            Var v = var(l);
            if (c.size() < occurrences[v].clauseSize) {
                occurrences[v].clauseSize = c.size();
                occurrences[v].count = 1;

                if (order_heap.inHeap(l.x)) order_heap.decrease(l.x);
                if (order_heap.inHeap((~l).x)) order_heap.decrease((~l).x);
            } else if (c.size() == occurrences[v].clauseSize) {
                occurrences[v].count++;
                if (order_heap.inHeap(l.x)) order_heap.decrease(l.x);
                if (order_heap.inHeap((~l).x)) order_heap.decrease((~l).x);
            }
        }
    #endif
    }

    inline void PropagationQueue::handleEventRestarted(const struct BCPRLStats& stats) {
        if (ENABLE_PRIORITY_BCP_RL) {
            current_bcpmode = bcprlManager.selectNextMode(current_bcpmode, stats);
            bcprlManager.clearScores(stats);
        }
    }

    inline void PropagationQueue::handleEventLearntClause(uint64_t lbd) {
        bcprlManager.handleEventLearntClause(lbd);
    }
}

#endif