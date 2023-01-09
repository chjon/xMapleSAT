/*******************************************************************************[AssignmentTrail.h]
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

#ifndef Minisat_AssignmentTrail_h
#define Minisat_AssignmentTrail_h

#include "core/SolverTypes.h"

namespace Minisat {
    // Forward declarations
    class Solver;

    /**
     * @brief This class keeps a record of the variable assignment order, the decision levels, and
     * the reasons for each of the assignments. 
     * 
     */
    class AssignmentTrail {
    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // HELPER TYPES

        /// @brief A wrapper to store the reason and level for a variable
        struct VarData {
            /// @brief The reason clause for the variable assignment
            CRef reason;

            /// @brief The decision level at which the variable was assigned
            int level;
        };

    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // MEMBER VARIABLES

        /// @brief The current variable assignments
        vec<lbool> assigns;

        /// @brief Assignment record; stores all assigments made in the order they were made.
        vec<Lit> trail;

        /// @brief Separator indices for different decision levels in 'trail'.
        vec<int> trail_lim;

        /// @brief Reason and level for each variable.
        vec<VarData> vardata;

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // SOLVER REFERENCES

        ClauseAllocator& ca;
        Solver& solver;

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////////
        // TEMPORARY VARIABLES

        /// @brief The total number of calls to LBD
        uint64_t lbd_calls;

        /// @brief Used to keep track of variables when computing LBD
        vec<uint64_t> lbd_seen;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CONSTRUCTORS

        /**
         * @brief Construct a new AssignmentTrail object
         * 
         * @param s Reference to main solver object
         */
        AssignmentTrail(Solver& s);

        /**
         * @brief Destroy the AssignmentTrail object
         * 
         */
        ~AssignmentTrail() = default;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // ACCESSORS

        /**
         * @brief Get the current number of variables
         * 
         * @return the current number of variables
         */
        int nVars (void) const;

        /**
         * @brief Get the truth assignment of a variable
         * 
         * @param x the variable whose truth assignment should be returned
         * @return l_True, l_False, or l_Undef, depending on the variable assignment
         */
        lbool value (Var x) const;

        /**
         * @brief Get the truth assignment of a literal
         * 
         * @param x the literal whose truth assignment should be returned
         * @return l_True, l_False, or l_Undef, depending on the underlying variable assignment
         */
        lbool value (Lit p) const;

        /**
         * @brief Get a view-only reference to the trail
         * 
         * @return a view-only reference to the trail
         */
        const vec<Lit>& getTrail(void) const;

        /**
         * @brief Index into the assignment trail
         * 
         * @param i the index for which to get the literal
         * @return the literal at the given index 
         */
        Lit operator[](int i) const;

        /**
         * @brief Get the first index of a decision level on the trail.
         * 
         * @param level the decision level for which to get the index.
         * @return The first index of the decision level on the trail.
         */
        int indexOfDecisionLevel(int level) const;

        /**
         * @brief Get the current decision level
         * 
         * @return The current decision level
         */
        int decisionLevel(void) const;

        /**
         * @brief Return an abstraction of sets of decision levels.
         * 
         * @param x The variable for which to get an abstract decision level.
         * @return A value representing the decision level of the variable 
         */
        uint32_t abstractLevel(Var x) const;

        /**
         * @brief Get the reason for a variable assignment. Assumes that x is assigned.
         * 
         * @param x The variable for which to get the reason.
         * @return The reason clause responsible for the variable assignment.
         */
        CRef reason(Var x) const;

        /**
         * @brief Get the decision level of a variable assignment. Assumes that x is assigned.
         * 
         * @param x The variable for which to get the decision level.
         * @return The decision level at which the variable was assigned.
         */
        int level(Var x) const;

        /**
         * @brief Get the current number of assigned literals.
         * 
         * @return the current number of assigned literals. 
         */
        int nAssigns(void) const;

        /**
         * @brief Get the current number of literals assigned at the root level.
         * 
         * @return the current number of literals assigned at the root level. 
         */
        int nRootAssigns(void) const;

        /**
         * @brief Determine whether a clause is a reason for some implication in the current state.
         * 
         * @param c the clause to check. Assumes that the clause variable ordering satisfies the
         * watcher invariant.
         * @return true if the clause is reason for some implication in the current state, false
         * otherwise.
         */
        bool locked(const Clause& c) const;

        /**
         * @brief Check whether a clause is satisfied under the current variable assignment
         * 
         * @param c the clause to check
         * @return true iff a literal in the clause is satisfied
         */
        bool satisfied(const Clause& c) const;

        double progressEstimate () const; // DELETE THIS ?? IT'S NOT VERY USEFUL ...

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATE MODIFICATION

        /**
         * @brief Set up internal data structures for a new variable
         * 
         * @param v the variable to register
         * @note Assumes that @code{v} has not already been registered
         * 
         * @return the ID of the new variable
         */
        Var newVar(void);

        /**
         * @brief Begins a new decision level.
         * 
         * @note Implemented by storing the index of the assignment trail where the new decision
         * level begins.
         */
        void newDecisionLevel(void);

        /**
         * @brief Assign a variable such that the given literal is true. 
         * @warning This might not add the literal to the propagation queue! Prefer to use the
         * PropagationQueue's @code{enqueue} method instead.
         * 
         * @param p The literal to assign. Assumes that the current value of the literal is
         * undefined.
         * @param from The reason for the literal assignment.
         */
        void assign(Lit p, CRef from = CRef_Undef);

        void simpleAssign(Lit p, CRef from = CRef_Undef);

        /**
         * @brief Backtrack until a certain decision level. Keeps all assignments at 'level' but
         * not beyond.
         * 
         * @param level The decision level to which to backtrack.
         */
        void cancelUntilLevel(int level);

        void cancelUntilTrailSize(int trailSize);

        /**
         * @brief Clean up when a clause is deleted - don't leave pointers to free'd memory!
         * 
         * @param c the clause that was deleted
         */
        void handleEventClauseDeleted(const Clause& c);

        /**
         * @brief Relocate all clauses
         * 
         * @param to the ClauseAllocator to relocate to
         */
        void relocAll(ClauseAllocator& to);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////////
        // UTILITY FUNCTIONS

        template<class V>
        int computeLBD(const V& c);

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////////
        // HELPER FUNCTIONS

        /**
         * @brief Assign a variable such that the given literal is true. 
         * @warning This might not add the literal to the propagation queue! Prefer to use the
         * PropagationQueue's @code{enqueue} method instead.
         * 
         * @tparam simple true to skip notifying event listeners
         * @param p The literal to assign. Assumes that the current value of the literal is
         * undefined.
         * @param from The reason for the literal assignment.
         */
        template <bool simple>
        void genericAssign(Lit p, CRef from);

        template<bool notifyListeners>
        void cancelUntil(int trailSize);
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE METHODS

    //////////////
    // ACCESSORS

    inline int AssignmentTrail::nVars() const {
        return assigns.size();
    }

    inline lbool AssignmentTrail::value(Var x) const {
        return assigns[x];
    }

    inline lbool AssignmentTrail::value(Lit p) const {
        return assigns[var(p)] ^ sign(p);
    }

    inline const vec<Lit>& AssignmentTrail::getTrail(void) const {
        return trail;
    }

    inline Lit AssignmentTrail::operator[](int i) const {
        return trail[i];
    }

    inline int AssignmentTrail::indexOfDecisionLevel(int l) const {
        return l == 0 ? 0 : trail_lim[l - 1];
    }

    inline int AssignmentTrail::decisionLevel() const {
        return trail_lim.size();
    }

    inline uint32_t AssignmentTrail::abstractLevel(Var x) const {
        return 1 << (level(x) & 31);
    }

    inline CRef AssignmentTrail::reason(Var x) const {
        return vardata[x].reason;
    }

    inline int AssignmentTrail::level(Var x) const {
        return vardata[x].level;
    }

    inline int AssignmentTrail::nAssigns() const {
        return trail.size();
    }

    inline int AssignmentTrail::nRootAssigns() const {
        return trail_lim.size() == 0 ? trail.size() : trail_lim[0];
    }

    inline bool AssignmentTrail::locked(const Clause& c) const {
        const int i = c.size() != 2 ? 0 : (value(c[0]) == l_True ? 0 : 1);
        return
            value(c[i]) == l_True &&
            reason(var(c[i])) != CRef_Undef &&
            ca.lea(reason(var(c[i]))) == &c;
    }

    inline bool AssignmentTrail::satisfied(const Clause& c) const {
        for (int i = 0; i < c.size(); i++)
            if (value(c[i]) == l_True) return true;
        return false;
    }

    ///////////////////////
    // STATE MODIFICATION

    inline Var AssignmentTrail::newVar(void) {
        const Var v = nVars();
        assigns.push(l_Undef);
        vardata.push(VarData{CRef_Undef, 0});
        trail  .capacity(v + 1);
        lbd_seen.push(0);
        return v;
    }

    inline void AssignmentTrail::newDecisionLevel(void) {
        trail_lim.push(trail.size());
    }

    inline void AssignmentTrail::handleEventClauseDeleted(const Clause& c) {
        if (locked(c)) {
            Lit implied = c.size() != 2 ? c[0] : (value(c[0]) == l_True ? c[0] : c[1]);
            vardata[var(implied)].reason = CRef_Undef;
        }
    }

    inline void AssignmentTrail::relocAll(ClauseAllocator& to) {
        for (int i = 0; i < trail.size(); i++) {
            Var v = var(trail[i]);
            if (reason(v) != CRef_Undef && (ca[reason(v)].reloced() || locked(ca[reason(v)])))
                ca.reloc(vardata[v].reason, to);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // UTILITY FUNCTIONS

    template <class V>
    int AssignmentTrail::computeLBD(const V& c) {
        int lbd = 0;
        lbd_calls++;
        for (int i = 0; i < c.size(); i++) {
            int l = level(var(c[i]));
            if (l != 0 && lbd_seen[l] != lbd_calls){
                lbd_seen[l] = lbd_calls;
                lbd++;
            }
        }

        return lbd;
    }
}

#endif