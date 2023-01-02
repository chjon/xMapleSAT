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

#include "core/SolverTypes.h"
#include "core/VariableDatabase.h"
#include "core/AssignmentTrail.h"
#include "core/PropagationQueue.h"
#include <mtl/Heap.h>
#include <map>

// Making some internal methods visible for testing
#ifdef TESTING
#define protected public
#endif

namespace Minisat {
    // Forward declarations
    class Solver;
    class ClauseDatabase;

    /**
     * @brief This class handles literal propagation.
     */
    class UnitPropagator {
    private:

        //////////////////////
        // HELPER FUNCTIONS //
        //////////////////////

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

        //////////////////////
        // MEMBER VARIABLES //
        //////////////////////

        // 'watches[lit]' is a list of constraints watching 'lit' (will go there if literal becomes true).
        OccLists<Lit, vec<Watcher>, WatcherDeleted> watches;

        //////////////////////////
        // RESOURCE CONSTRAINTS //
        //////////////////////////

        int64_t propagation_budget; // -1 means no budget.

    public:
        ////////////////
        // STATISTICS //
        ////////////////

        uint64_t propagations  ; // Total number of propagations performed by @code{propagate}

    protected:
        ///////////////////////
        // SOLVER REFERENCES //
        ///////////////////////

        VariableDatabase& variableDatabase;
        AssignmentTrail& assignmentTrail;
        PropagationQueue& propagationQueue;
        ClauseAllocator& ca;
        ClauseDatabase& clauseDatabase;
        Solver& solver;

    public:
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
         * @param c The clause to attach.
         * @param cr The CRef of the clause to attach. Must match @code{c}.
         * 
         * @note c and cr are provided separately to enable compiler optimization at the caller.
         */
        void attachClause(const Clause& c, CRef cr);

        /**
         * @brief Attach a clause to watcher lists.
         * 
         * @param cr The CRef of the clause to attach.
         */
        void attachClause(CRef cr);

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

        /**
         * @brief Get list of non-binary clause watchers for a given literal
         * 
         * @param l the literal whose watchers should be returned
         * @return non-binary clause watchers for @code{l}
         */
        const vec<Watcher>& getWatchers(Lit l) const;

        /**
         * @brief Get list of binary clause watchers for a given literal
         * 
         * @param l the literal whose watchers should be returned
         * @return binary clause watchers for @code{l}
         */
        // const vec<Watcher>& getBinaryWatchers(Lit l) const;

        /**
         * @brief Move undefined literal to index 0, ensuring that watcher invariants are satisfied
         * 
         * @param cr The CRef of the asserting clause
         * @param i_undef The index of the undefined literal in the clause
         * @param i_max The index of the literal in the clause with the highest decision level
         */
        void enforceWatcherInvariant(CRef cr, int i_undef, int i_max);

        //////////////////////////
        // RESOURCE CONSTRAINTS //
        //////////////////////////

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

        /**
         * @brief Check whether the solver can still propagate within the budget.
         * 
         * @return false if the solver has exceeded the budget, true otherwise
         */
        bool withinBudget() const;
    };

    inline void UnitPropagator::newVar(Var v) {
        watches.init(mkLit(v, false));
        watches.init(mkLit(v, true ));
    }

    inline void UnitPropagator::attachClause(const Clause& c, CRef cr) {
        assert(c.size() > 1);
        OccLists<Lit, vec<Watcher>, WatcherDeleted>& ws = watches;
        ws[~c[0]].push(Watcher(cr, c[1]));
        ws[~c[1]].push(Watcher(cr, c[0]));
    }
    inline void UnitPropagator::attachClause(CRef cr) { attachClause(ca[cr], cr); }

    inline void UnitPropagator::detachClause(const Clause& c, CRef cr, bool strict) {
        assert(c.size() > 1);
        OccLists<Lit, vec<Watcher>, WatcherDeleted>& ws = watches;
        if (strict) {
            remove(ws[~c[0]], Watcher(cr, c[1]));
            remove(ws[~c[1]], Watcher(cr, c[0]));
        } else {
            // Lazy detaching: (NOTE! Must clean all watcher lists before garbage collecting this clause)
            ws.smudge(~c[0]);
            ws.smudge(~c[1]);
        }
    }
    inline void UnitPropagator::detachClause(CRef cr, bool strict) { detachClause(ca[cr], cr, strict); }

    inline void UnitPropagator::relocWatchers(vec<Watcher>& ws, ClauseAllocator& to) {
        for (int i = 0; i < ws.size(); i++) ca.reloc(ws[i].cref, to);
    }

    inline const vec<Watcher>& UnitPropagator::getWatchers(Lit l) const { return watches    [l]; }

    inline void UnitPropagator::setPropBudget(int64_t x) { propagation_budget = propagations + x; }
    inline void UnitPropagator::budgetOff() { propagation_budget = -1; }
    inline bool UnitPropagator::withinBudget() const { return propagation_budget < 0 || propagations < (uint64_t) propagation_budget; }

} // namespace Minisat

#endif