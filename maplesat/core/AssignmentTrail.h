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

#ifndef Minisat_AssignmentTrail_h
#define Minisat_AssignmentTrail_h

#include "core/SolverTypes.h"
#include "core/VariableDatabase.h"

namespace Minisat {
    class Solver;
    class PropagationQueue;

    /**
     * @brief This class keeps a record of the variable assignment order, the decision levels, and
     * the reasons for each of the assignments. 
     * 
     */
    class AssignmentTrail {
    protected:
        /**
         * @brief A wrapper to store the reason and level for a variable
         * 
         */
        struct VarData {
            CRef reason; // The reason clause for the variable assignment
            int level;   // The decision level at which the variable was assigned
        };

        /**
         * @brief Constructor for the VarData object.
         * 
         * @param cr The reason clause for the variable assignment
         * @param l The decision level at which the variable was assigned
         * @return A VarData object containing the arguments.
         */
        static inline VarData mkVarData(CRef cr, int l) { VarData d = {cr, l}; return d; }

        // Assignment record; stores all assigments made in the order they were made.
        vec<Lit> trail;

        // Separator indices for different decision levels in 'trail'.
        vec<int> trail_lim;

        // Reason and level for each variable.
        vec<VarData> vardata;

    public:
        //////////////////
        // CONSTRUCTORS //
        //////////////////

        AssignmentTrail(Solver* s);
        ~AssignmentTrail() = default;

        ////////////////
        // PUBLIC API //
        ////////////////

        /**
         * @brief Set up internal data structures for a new variable
         * 
         * @param v the variable to register
         * @note Assumes that @code{v} has not already been registered
         */
        void newVar(Var v);

        /**
         * @brief Relocate all clauses
         * 
         * @param to the ClauseAllocator to relocate to
         */
        void relocAll(ClauseAllocator& to);

        /**
         * @brief Begins a new decision level.
         * 
         * @note Implemented by storing the index of the assignment trail where the new decision level begins
         */
        void newDecisionLevel(void);

        /**
         * @brief Assign a variable such that the given literal is true. 
         * @warning This might not add the literal to the propagation queue! Prefer to use the
         * PropagationQueue's @code{enqueue} method instead.
         * 
         * @param p The literal to assign. Assumes that the current value of the literal is undefined.
         * @param from The reason for the literal assignment.
         */
        void assign(Lit p, CRef from = CRef_Undef);

        /**
         * @brief Backtrack until a certain decision level.
         * 
         * @param level The decision level to which to backtrack.
         */
        void cancelUntil(int level);

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
         * @brief Get the reason for a variable assignment. Assumes that the variable is assigned.
         * 
         * @param x The variable for which to get the reason.
         * @return The reason clause responsible for the variable assignment.
         */
        CRef reason(Var x) const;

        /**
         * @brief Get the decision level of a variable assignment. Assumes that the variable is assigned.
         * 
         * @param x The variable for which to get the decision level.
         * @return The decision level at which the variable was assigned.
         */
        int level(Var x) const;

        /**
         * @brief Determine whether a clause is a reason for some implication in the current state.
         * 
         * @param c the clause to check. Assumes that the clause variable ordering satisfies the watcher invariant.
         * @return true if the clause is reason for some implication in the current state, false otherwise.
         */
        bool locked (const Clause& c) const;

        void handleEventClauseDeleted(const Clause& c);

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
         * @brief Index into the assignment trail
         * 
         * @param i the index for which to get the literal
         * @return the literal at the given index 
         */
        Lit operator[](int i) const;

        int indexOfDecisionLevel(int l) const;

        double progressEstimate () const; // DELETE THIS ?? IT'S NOT VERY USEFUL ...

    private:
        VariableDatabase& variableDatabase;
        ClauseAllocator& ca;
        Solver* solver;

        friend PropagationQueue;
    };

    inline void AssignmentTrail::newVar(Var v) {
        vardata.push(mkVarData(CRef_Undef, 0));
        trail  .capacity(v + 1);
    }

    inline void     AssignmentTrail::newDecisionLevel(void)             { trail_lim.push(trail.size()); }
    
    inline int      AssignmentTrail::decisionLevel   ()      const      { return trail_lim.size(); }
    inline uint32_t AssignmentTrail::abstractLevel   (Var x) const      { return 1 << (level(x) & 31); }
    inline CRef     AssignmentTrail::reason          (Var x) const      { return vardata[x].reason; }
    inline int      AssignmentTrail::level           (Var x) const      { return vardata[x].level; }
    
    inline bool AssignmentTrail::locked (const Clause& c) const { return variableDatabase.value(c[0]) == l_True && reason(var(c[0])) != CRef_Undef && ca.lea(reason(var(c[0]))) == &c; }
    inline void AssignmentTrail::handleEventClauseDeleted(const Clause& c) { if (locked(c)) vardata[var(c[0])].reason = CRef_Undef; }

    inline int AssignmentTrail::nAssigns  () const { return trail.size(); }
    inline int AssignmentTrail::nRootAssigns() const { return trail_lim.size() == 0 ? trail.size() : trail_lim[0]; }

    inline Lit AssignmentTrail::operator[] (int i) const { return trail[i]; }
    inline int AssignmentTrail::indexOfDecisionLevel(int l) const { return l == 0 ? 0 : trail_lim[l - 1]; }
}

#endif