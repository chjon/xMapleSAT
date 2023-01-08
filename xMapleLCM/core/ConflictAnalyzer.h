/******************************************************************************[ConflictAnalyzer.h]
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

#ifndef Minisat_ConflictAnalyzer_h
#define Minisat_ConflictAnalyzer_h

#include "core/SolverTypes.h"

namespace Minisat {
    // Forward declarations
    class Solver;
    class AssignmentTrail;
    class BranchingHeuristicManager;

    /**
     * @brief This class is responsible for analyzing conflict graphs to generate and simplify
     * learnt clauses.
     * 
     */
    class ConflictAnalyzer {
    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // HELPER TYPES

        enum ConflictClauseMinimizationMode: int {
            NONE  = 0,
            BASIC = 1,
            DEEP  = 2,
        };

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // SOLVER REFERENCES

        AssignmentTrail& assignmentTrail;
        BranchingHeuristicManager& branchingHeuristicManager;
        ClauseAllocator& ca;
        Solver& solver;

    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // PARAMETERS
        
        /// @brief Controls conflict clause minimization
        ConflictClauseMinimizationMode ccmin_mode;

        ///////////////////////////////////////////////////////////////////////////////////////////
        // TEMPORARY VARIABLES
        //
        // these variables are allocated here to avoid repeated allocation costs

        /// @brief A learnt clause
        vec<Lit> firstUIPClause;

        /// @brief Work stack for @code{litRedundant}: holds list of reason clauses to be examined
        vec<CRef> workStack;

        /// @brief Marks whether a variable has already been seen
        vec<bool> seen;

        /// @brief Mostly for efficient LBD computation. 'seen2[i]' will indicate if decision level or variable 'i' has been seen.
        vec<uint64_t> seen2;
        
        /// @brief Simple counter for marking purpose with 'seen2'.
        uint64_t counter;

        /// @brief Stores variables whose values have been set in @code{seen} and need to be
        /// cleared. Currently only used by @code{litRedundant}
        vec<Var> toClear;

        // CollectFirstUIP data structures

        vec<int> pathCs;

        vec<double> var_iLevel_tmp;

        vec<Lit> involved_lits;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATISTICS

        /// @brief The total number of literals in from first UIP learnt clauses
        uint64_t max_literals;

        /// @brief The total number of literals in learnt clauses after clause simplification and
        /// minimization
        uint64_t tot_literals;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CONSTRUCTORS

        /**
         * @brief Construct a new ConflictAnalyzer object
         * 
         * @param s Reference to main solver object
         */
        ConflictAnalyzer(Solver& s);

        /**
         * @brief Destroy the ConflictAnalyzer object
         * 
         */
        ~ConflictAnalyzer() = default;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATE MODIFICATION

        /**
         * @brief Set up internal data structures for a new variable
         * 
         * @param v the variable to register
         */
        void newVar(Var v);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // PUBLIC API

        /**
         * @brief Analyze the current conflict graph and generate a learnt clause and its
         * associated backtrack level.
         * 
         * @param confl The clause responsible for the current conflict
         * @param out_learnt Output: the learnt clause
         * @param out_btlevel Output: the backtrack level for the learnt clause
         * @param out_lbd Output: the LBD of the learnt clause
         * 
         * @pre 'out_learnt' is assumed to be cleared
         * @pre Current decision level must be greater than root level
         * 
         * @post 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
         * @post If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of
         * the rest of the literals. There may be others from the same level.
         */
        void analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel, int& out_lbd);

        void simpleAnalyze(CRef confl, vec<Lit>& out_learnt, vec<CRef>& reason_clause, bool True_confl, int trailRecord);

        /**
         * @brief Specialized analysis procedure to express the final conflict in terms of
         * assumptions. Calculates the (possibly empty) set of assumptions that led to the
         * assignment of 'p'.
         * 
         * @param p the conflicting literal
         * @param out_conflict the set of assumptions leading to the assignment of @code{p}
         */
        void analyzeFinal(Lit p, vec<Lit>& out_conflict);

        bool collectFirstUIP(CRef confl);

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // HELPER FUNCTIONS

        /**
         * @brief Check whether a literal is redundant and can be removed.
         * 
         * @param p the literal to check for redundancy
         * @param abstract_levels an abstraction of decision levels, used to abort early if the
         * algorithm is visiting literals at levels that cannot be removed later.
         * @return true if p is redundant and can be removed, false otherwise
         * 
         * @note this is a helper method for @code{simplifyClauseDeep}
         * @pre the @code{workStack} vector is empty
         * @post the @code{workStack} vector is empty
         * 
         */
        bool litRedundant(Lit p, uint32_t abstract_levels);

        /**
         * @brief Simplify a learnt clause
         * 
         * @param simplified the output simplified clause.
         * @param toSimplify the clause to simplify.
         */
        void simplifyClauseDeep(vec<Lit>& simplified, const vec<Lit>& toSimplify);

        /**
         * @brief Check whether a reason clause is subsumed by a set of variables
         * 
         * @param c the reason clause for a variable in the learnt clause
         * @return true if the reason clause is subsumed, false otherwise
         * 
         * @note This is a helper method for @code{simplifyClauseBasic}
         * @note Precondition: @code{seen} is true for the variables in the learnt clause
         */
        bool reasonSubsumed(const Clause& c);

        /**
         * @brief Simplify a learnt clause
         * 
         * @param simplified the output simplified clause.
         * @param learntClause the learnt clause to modify.
         */
        void simplifyClauseBasic(vec<Lit>& simplified, const vec<Lit>& toSimplify);

        /**
         * @brief Learn a clause according to the first UIP learning scheme.
         * @param confl the conflict clause
         * @param learntClause the output learnt clause. Assumed to be empty initially.
         */
        void getFirstUIPClause(CRef confl, vec<Lit>& learntClause);

        /**
         * @brief Simplify a learnt clause
         * 
         * @param simplified the output simplified clause.
         * @param toSimplify the learnt clause to modify.
         */
        void simplifyClause(vec<Lit>& simplified, const vec<Lit>& toSimplify);

        /**
         * @brief Further learnt clause minimization by binary resolution.
         * 
         * @param out_learnt the learnt clause to modify.
         * @return true iff the learnt clause was modified.
         */
        bool binResMinimize(vec<Lit>& out_learnt);

        /**
         * @brief Enforce the watcher invariant for a learnt clause. Assumes that the variable with
         * the highest decision level is at index 0. Moves the variable with the next-highest
         * decision level to index 1.
         * 
         * @param learntClause the learnt clause to modify.
         */
        void enforceWatcherInvariant(vec<Lit>& learntClause);
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE METHODS

    inline void ConflictAnalyzer::newVar(Var v) {
        seen.push(0);
        seen2.push(0);
        pathCs.push(0);
        var_iLevel_tmp.push(0);
    }
}

#endif