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
#define BCP_PRIORITY_OUT_OF_ORDER 2 // Out-of-order propagation

// Select prioritization mode
#ifndef BCP_PRIORITY_MODE
    #define BCP_PRIORITY_MODE BCP_PRIORITY_IMMEDIATE
#endif

#include "core/AssignmentTrail.h"
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
         * @brief Propagation modes for priority BCP
         * 
         */
        enum PropagationMode: int {
            IMMEDIATE    = BCP_PRIORITY_IMMEDIATE,
            DELAYED      = BCP_PRIORITY_DELAYED,
            OUT_OF_ORDER = BCP_PRIORITY_OUT_OF_ORDER,
        };

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
            LitOrderLt(const vec<T>&  act) : activity(&act) { }
        };

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // SOLVER REFERENCES

        AssignmentTrail& assignmentTrail;

    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // MEMBER VARIABLES

        /// @brief The propagation mode to use
        PropagationMode propagationMode;

        /// @brief The head of the propagation queue as an index into the assignment trail
        int qhead;

        /// @brief A view-only reference to the assignment trail
        const vec<Lit>& queue;

        /// @brief The priority queue for selecting variables to propagate during BCP
        Heap< LitOrderLt<double> > order_heap;

        /// @brief The polarities of queued variables
        vec<lbool> soft_assigns;

        /// @brief The reason clauses for queuing variables for BCP
        vec<CRef> reasons;

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

        /**
         * @brief Add a literal to the propagation queue without notifying event listeners
         * 
         * @param p the literal to add to the queue
         * @param from the reason for the literal
         * @return false if the negated version of the literal is already in the queue, true otherwise
         */
        bool simpleEnqueue(Lit p, CRef from = CRef_Undef);

        /**
         * @brief Add a set of literals to the propagation queue
         * 
         * @param trail the current assignment trail
         * @param levelHead the start index from which to add literals to the queue
         * 
         * @note literals are added from @code{trail[levelHead]} up to @code{trail.last()}
         */
        void batchEnqueue(const vec<Lit>& trail, int levelHead);

        /**
         * @brief Get the next literal for propagation
         * 
         * @tparam simple: true to skip updating stats and notifying event listeners
         * @return the literal to propagate
         */
        template <bool simple>
        Lit getNext();

        /**
         * @brief Clear the propagation queue
         * 
         */
        void clear();

        /**
         * @brief Set the prioritization mode and activity metric to use for priority BCP
         * 
         * @param activity the activities of each variable
         */
        void updatePrioritizationMode(
            bool distanceBranching,
            const vec<double>& activity
        );
    
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
        template <bool simple>
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
    }

    inline bool PropagationQueue::enqueue(Lit p, CRef from) {
        return genericEnqueue<false>(p, from);
    }

    inline bool PropagationQueue::simpleEnqueue(Lit p, CRef from) {
        return genericEnqueue<true>(p, from);
    }

    inline void PropagationQueue::batchEnqueue(const vec<Lit>& trail, int levelHead) {
        switch (propagationMode) {
            case PropagationMode::IMMEDIATE:
            case PropagationMode::DELAYED: {
                qhead = levelHead;
            } break;

            case PropagationMode::OUT_OF_ORDER: {
                const Lit* i; const Lit* end;
                i = static_cast<const Lit*>(trail);
                end = i + trail.size();
                while (i != end) {
                    order_heap.insert((i++)->x);
                }
            } break;

            default: assert(false); // Propagation mode must be one of the above
        }
    }

    template <bool simple>
    inline Lit PropagationQueue::getNext() {
        switch (propagationMode) {
            case PropagationMode::IMMEDIATE: {
                return qhead == queue.size() ? lit_Undef : queue[qhead++];
            }

            case PropagationMode::DELAYED: {
                // Propagate along the trail first
                if (qhead < queue.size()) return queue[qhead++];

                // Propagate according to priority queue afterward
                qhead++;
                if (!order_heap.empty()) {
                    Lit p = Lit{order_heap.removeMin()};

                    // Note: Variable is assigned here because it is not assigned in @code{enqueue()}
                    if (simple) assignmentTrail.simpleAssign(p, reasons[var(p)]);
                    else        assignmentTrail.assign(p, reasons[var(p)]);
                    soft_assigns[var(p)] = l_Undef;
                    return p;
                }
                return lit_Undef;
            }

            case PropagationMode::OUT_OF_ORDER: {
                // Always propagate according to priority queue
                return order_heap.empty() ? lit_Undef : Lit{order_heap.removeMin()};
            }

            default: {
                assert(false);
                return lit_Undef; // Propagation mode must be one of the above
            }
        }
    }

    inline void PropagationQueue::clear() {
        switch (propagationMode) {
            case PropagationMode::IMMEDIATE: {
                qhead = queue.size();
            } break;

            case PropagationMode::DELAYED: {
                qhead = queue.size();
                for (int k = 0; k < order_heap.size(); k++)
                    soft_assigns[order_heap[k] >> 1] = l_Undef;
                order_heap.clear();
            } break;

            case PropagationMode::OUT_OF_ORDER: {
                order_heap.clear();
            } break;

            default: assert(false); // Propagation mode must be one of the above
        }
    }

    inline void PropagationQueue::updatePrioritizationMode(
        bool distanceBranching,
        const vec<double>& activity
    ) {
        PropagationMode prevPropagationMode = propagationMode;
        propagationMode = distanceBranching
            ? PropagationMode::IMMEDIATE
            : static_cast<PropagationMode>(BCP_PRIORITY_MODE);

        // Check whether prioritization mode changed
        if (propagationMode != prevPropagationMode) {
            batchEnqueue(queue, 0);
        }

        // Set comparator
        switch (propagationMode) {
            case PropagationMode::DELAYED:
            case PropagationMode::OUT_OF_ORDER: {
                order_heap.setComp(LitOrderLt<double>(activity));
            } break;

            case PropagationMode::IMMEDIATE:
            default: break;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // HELPER FUNCTIONS

    template<bool simple>
    inline bool PropagationQueue::genericEnqueue(Lit p, CRef from) {
        switch (propagationMode) {
            case PropagationMode::IMMEDIATE: {
                if (simple) assignmentTrail.simpleAssign(p, from);
                else        assignmentTrail.assign(p, from);
            } break;
            case PropagationMode::DELAYED: {
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
            case PropagationMode::OUT_OF_ORDER: {
                order_heap.insert(p.x);
                if (simple) assignmentTrail.simpleAssign(p, from);
                else        assignmentTrail.assign(p, from);
            } break;

            default: break;
        }

        return true;
    }
}

#endif