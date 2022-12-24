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

#include "core/SolverTypes.h"
#include "core/AssignmentTrail.h"
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
        ////////////////////
        // HELPER STRUCTS //
        ////////////////////

        // Comparator for BCP priority queue
        template<class T>
        struct LitOrderLt {
            const vec<T>&  activity;
            bool operator () (Var x, Var y) const {
                x >>= 1; y >>= 1;
                return activity[x] > activity[y];
            }
            LitOrderLt(const vec<T>&  act) : activity(act) { }
        };

        //////////////////////
        // MEMBER VARIABLES //
        //////////////////////

    #if BCP_PRIORITY_MODE == BCP_PRIORITY_IMMEDIATE
        int qhead;
        vec<Lit>& queue;

    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
        int qhead;
        vec<Lit>& queue;
        Heap< LitOrderLt<double> > order_heap;
        vec<lbool> soft_assigns;
        vec<CRef>  reasons;

    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
        Heap< LitOrderLt<double> > order_heap;

    #endif
        VariableDatabase& variableDatabase;
        AssignmentTrail& assignmentTrail;

        //////////////////////
        // HELPER FUNCTIONS //
        //////////////////////

    public:
        //////////////////
        // CONSTRUCTORS //
        //////////////////

        PropagationQueue(Solver* s);
        ~PropagationQueue() = default;

        ////////////////
        // PUBLIC API //
        ////////////////

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
         * @return the literal to propagate
         */
        Lit getNext();

        /**
         * @brief Clear the propagation queue
         * 
         */
        void clear();
    };

    // Explicitly instantiate required templates
    template class PropagationQueue::LitOrderLt<double>;

    ////////////////////////////////////////
    // IMPLEMENTATION OF INLINE FUNCTIONS //
    ////////////////////////////////////////

    inline void PropagationQueue::newVar(Var v) {
    #if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
        soft_assigns.push(l_Undef);
        reasons.push(CRef_Undef);
    #endif
    }

    inline bool PropagationQueue::enqueue(Lit p, CRef from) {
    #if BCP_PRIORITY_MODE == BCP_PRIORITY_IMMEDIATE
        assignmentTrail.uncheckedEnqueue(p, from);

    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
        if (soft_assigns[var(p)] ^ sign(p) == l_False) {
            // Ensure conflicting literal is on the trail
            if (variableDatabase.value(p) == l_Undef)
                assignmentTrail.uncheckedEnqueue(~p, reasons[var(p)]);

            return false;
        } else if (soft_assigns[var(p)] ^ sign(p) == l_Undef) {
            // Note: actual variable assignment is delayed until the literal is popped by @code{getNext()}
            order_heap.insert(p.x);
            soft_assigns[var(p)] = lbool(!sign(p));
            reasons[var(p)] = from;
        }

    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
        order_heap.insert(p.x);
        assignmentTrail.uncheckedEnqueue(p, from);

    #endif

        return true;
    }

    inline void PropagationQueue::batchEnqueue(vec<Lit>& trail, int levelHead) {
    #if BCP_PRIORITY_MODE == BCP_PRIORITY_IMMEDIATE
        qhead = levelHead;

    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
        qhead = levelHead;
        
    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
        Lit* i, end;
        i = static_cast<Lit*>(trail);
        end = i + trail.size();
        while (i != end) {
            order_heap.insert((i++)->x);
        }
    #endif
    }

    inline Lit PropagationQueue::getNext() {
    #if BCP_PRIORITY_MODE == BCP_PRIORITY_IMMEDIATE
        return qhead == queue.size() ? lit_Undef : queue[qhead++];
    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
        // Propagate along the trail first
        if (qhead < queue.size()) return queue[qhead++];

        // Propagate according to priority queue afterward
        qhead++;
        if (!order_heap.empty()) {
            Lit p = Lit(order_heap.removeMin());

            // Note: Variable is assigned here because it is not assigned in @code{enqueue()}
            assignmentTrail.uncheckedEnqueue(p, from);
            return p;
        }
        return lit_Undef;

    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
        // Always propagate according to priority queue
        return order_heap.empty() ? lit_Undef : Lit(order_heap.removeMin());
    #endif
    }

    inline void PropagationQueue::clear() {
    #if BCP_PRIORITY_MODE == BCP_PRIORITY_IMMEDIATE
        qhead = queue.size();

    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
        qhead = queue.size();
        for (int k = 0; k < order_heap.size(); k++)
            soft_assigns[order_heap[k] >> 1] = l_Undef;
        order_heap.clear();

    #elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
        order_heap.clear();

    #endif
    }
}

#endif