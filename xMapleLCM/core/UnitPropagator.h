/********************************************************************************[UnitPropagator.h]
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

#ifndef Minisat_UnitPropagator_h
#define Minisat_UnitPropagator_h

#ifdef TESTING
#define protected public
#define private public
#endif

#include <map>
#include "core/SolverTypes.h"
#include "mtl/Heap.h"

namespace Minisat {
    // Forward declarations
    class Solver;
    class PropagationQueue;
    class AssignmentTrail;

    /**
     * @brief This class handles literal propagation.
     */
    class UnitPropagator {
    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // SOLVER REFERENCES

        PropagationQueue& propagationQueue;
        AssignmentTrail& assignmentTrail;
        ClauseAllocator& ca;
        Solver& solver;

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // MEMBER VARIABLES

        /// @brief is a list of binary constraints watching 'lit' (will go there if literal becomes
        /// true).
        OccLists<Lit, vec<Watcher>, WatcherDeleted> watches_bin;

        /// @brief is a list of non-binary constraints watching 'lit' (will go there if literal
        /// becomes true).
        OccLists<Lit, vec<Watcher>, WatcherDeleted> watches;

    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // PARAMETERS

        /// @brief The number of allowed propagations. -1 means no budget.
        int64_t propagation_budget;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATISTICS

        /// @brief Total number of propagations performed by @code{propagate}
        uint64_t propagations;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CONSTRUCTORS

        /**
         * @brief Construct a new UnitPropagator object
         * 
         * @param s Reference to main solver object
         */
        UnitPropagator(Solver& s);

        /**
         * @brief Destroy the UnitPropagator object
         */
        ~UnitPropagator() = default;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // ACCESSORS
        
        /**
         * @brief Get list of binary clause watchers for a given literal
         * 
         * @param l the literal whose binary watchers should be returned
         * @return clause watchers for @code{l}
         */
        const vec<Watcher>& getBinWatchers(Lit l) const;

        /**
         * @brief Get list of non-binary clause watchers for a given literal
         * 
         * @param l the literal whose non-binary watchers should be returned
         * @return clause watchers for @code{l}
         */
        const vec<Watcher>& getWatchers(Lit l) const;

        /**
         * @brief Check whether the solver can still propagate within the budget.
         * 
         * @return false if the solver has exceeded the budget, true otherwise
         */
        bool withinBudget() const;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // UTILITY FUNCTIONS

        /**
         * @brief Move undefined literal to index 0, ensuring that watcher invariants are satisfied
         * 
         * @param cr The CRef of the asserting clause
         * @param i_undef The index of the undefined literal in the clause
         * @param i_max The index of the literal in the clause with the highest decision level
         */
        void enforceWatcherInvariant(CRef cr, int i_undef, int i_max);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // PARAMETER MODIFICATION

        /**
         * @brief Set the propagation budget.
         * 
         * @param x The number of times left to propagate.
         */
        void setPropBudget(int64_t x);

        /**
         * @brief Set the solver to ignore the propagation budget. 
         * 
         */
        void budgetOff();

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATE MODIFICATION

        /**
         * @brief Register watchers for a new variable
         * 
         * @param v the variable to register
         * @note Assumes that @code{v} has not already been registered
         */
        void newVar(Var v);

        /**
         * @brief Attach a clause to watcher lists.
         * 
         * @param cr The CRef of the clause to attach.
         */
        void attachClause(CRef cr);

        /**
         * @brief Detach a clause from watcher lists.
         * 
         * @param cr The CRef of the clause to detach.
         * @param strict False to use lazy detaching, true otherwise
         */
        void detachClause(CRef cr, bool strict = false);

        /**
         * @brief Relocate CRefs to new ClauseAllocator
         * 
         * @param to The ClauseAllocator into which to reloc 
         */
        void relocAll(ClauseAllocator& to);

        /**
         * @brief Propagate all enqueued facts.
         * 
         * @return The conflicting clause if a conflict arises, otherwise CRef_Undef. 
         */
        CRef propagate();

        CRef simplePropagate();

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // HELPER FUNCTIONS

        /**
         * @brief Attach a clause to watcher lists.
         * 
         * @param c The clause to attach.
         * @param cr The CRef of the clause to attach. Must match @code{c}.
         * 
         * @note c and cr are provided separately to enable compiler optimization at the caller.
         */
        void attachClause(const Clause& c, CRef cr);

        /**
         * @brief Detach a clause from watcher lists.
         * 
         * @param c The clause to detach.
         * @param cr The CRef of the clause to detach. Must match @code{c}.
         * @param strict False to use lazy detaching, true otherwise
         * 
         * @note c and cr are provided separately to enable compiler optimization at the caller.
         */
        void detachClause(const Clause& c, CRef cr, bool strict = false);

        /**
         * @brief Perform all binary-clause propagations for a single literal
         * 
         * @param p the literal to propagate
         * @return The conflicting clause if a conflict arises, otherwise CRef_Undef.
         */
        CRef propagateSingleBinary(Lit p);

        /**
         * @brief Perform all non-binary-clause propagations for a single literal
         * 
         * @param p the literal to propagate
         * @return The conflicting clause if a conflict arises, otherwise CRef_Undef.
         */
        CRef propagateSingleNonBinary(Lit p);

        /**
         * @brief Perform all propagations for a single literal
         * 
         * @param p the literal to propagate
         * @return The conflicting clause if a conflict arises, otherwise CRef_Undef.
         */
        CRef propagateSingle(Lit p);

        /**
         * @brief Relocate watcher CRefs to new ClauseAllocator
         * 
         * @param ws The watchers for which to relocate CRefs
         * @param to The ClauseAllocator into which to reloc 
         */
        void relocWatchers(vec<Watcher>& ws, ClauseAllocator& to);
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE METHODS

    //////////////
    // ACCESSORS

    inline const vec<Watcher>& UnitPropagator::getBinWatchers(Lit l) const {
        return watches_bin[l];
    }

    inline const vec<Watcher>& UnitPropagator::getWatchers(Lit l) const {
        return watches[l];
    }

    inline bool UnitPropagator::withinBudget() const {
        return
            propagation_budget < 0 ||
            propagations < (uint64_t) propagation_budget;
    }

    ///////////////////////////
    // PARAMETER MODIFICATION

    inline void UnitPropagator::setPropBudget(int64_t x) {
        propagation_budget = propagations + x;
    }

    inline void UnitPropagator::budgetOff() {
        propagation_budget = -1;
    }

    ///////////////////////
    // STATE MODIFICATION

    inline void UnitPropagator::newVar(Var v) {
        watches_bin.init(mkLit(v, false));
        watches_bin.init(mkLit(v, true ));
        watches    .init(mkLit(v, false));
        watches    .init(mkLit(v, true ));
    }

    inline void UnitPropagator::attachClause(CRef cr) {
        attachClause(ca[cr], cr);
    }

    inline void UnitPropagator::detachClause(CRef cr, bool strict) {
        detachClause(ca[cr], cr, strict);
    }

    /////////////////////
    // HELPER FUNCTIONS

    inline void UnitPropagator::attachClause(const Clause& c, CRef cr) {
        assert(c.size() > 1);
        OccLists<Lit, vec<Watcher>, WatcherDeleted>& ws = (c.size() == 2) ? watches_bin : watches;
        ws[~c[0]].push(Watcher(cr, c[1]));
        ws[~c[1]].push(Watcher(cr, c[0]));
    }

    inline void UnitPropagator::detachClause(const Clause& c, CRef cr, bool strict) {
        assert(c.size() > 1);
        OccLists<Lit, vec<Watcher>, WatcherDeleted>& ws = (c.size() == 2) ? watches_bin : watches;
        if (strict) {
            remove(ws[~c[0]], Watcher(cr, c[1]));
            remove(ws[~c[1]], Watcher(cr, c[0]));
        } else {
            // Lazy detaching:
            // (NOTE! Must clean all watcher lists before garbage collecting this clause)
            ws.smudge(~c[0]);
            ws.smudge(~c[1]);
        }
    }

    inline void UnitPropagator::relocWatchers(vec<Watcher>& ws, ClauseAllocator& to) {
        for (int i = 0; i < ws.size(); i++) ca.reloc(ws[i].cref, to);
    }
} // namespace Minisat

#endif