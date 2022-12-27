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

#ifndef Minisat_ConflictAnalyzer_h
#define Minisat_ConflictAnalyzer_h

#include "core/SolverTypes.h"
#include "core/AssignmentTrail.h"
#include "core/BranchingHeuristicManager.h"

namespace Minisat {
    // Forward declarations
    class Solver;

    class ConflictAnalyzer {
    protected:
        ////////////////
        // PARAMETERS //
        ////////////////
        enum ConflictClauseMinimizationMode: int {
            NONE  = 0,
            BASIC = 1,
            DEEP  = 2,
        };
        
        // Controls conflict clause minimization
        ConflictClauseMinimizationMode ccmin_mode;

        //////////////////////
        // MEMBER VARIABLES //
        //////////////////////

        vec<Lit> firstUIPClause;

        // Work stack for @code{litRedundant}: holds list of reason clauses to be examined
        vec<CRef> workStack;

        vec<Var> toClear; // Temporary variable -- only used by @code{litRedundant}

    public:
        ////////////////
        // STATISTICS //
        ////////////////

        uint64_t max_literals;
        uint64_t tot_literals;

    protected:
        ///////////////////////
        // SOLVER REFERENCES //
        ///////////////////////

        AssignmentTrail& assignmentTrail;
        BranchingHeuristicManager& branchingHeuristicManager;
        ClauseAllocator& ca;
        Solver* solver;

        //////////////////////
        // HELPER FUNCTIONS //
        //////////////////////

        /**
         * @brief Check whether a literal is redundant and can be removed.
         * 
         * @param p the literal to check for redundancy
         * @param abstract_levels an abstraction of decision levels, used to abort early if the algorithm
         * is visiting literals at levels that cannot be removed later.
         * @return true if p is redundant and can be removed, false otherwise
         * 
         * @note this is a helper method for @code{simplifyClauseDeep}
         */
        bool litRedundant(Lit p, uint32_t abstract_levels);

        /**
         * @brief Simplify a learnt clause
         * 
         * @param learntClause the learnt clause to modify.
         */
        void simplifyClauseDeep(vec<Lit>& learntClause);

        /**
         * @brief Check whether a reason clause is subsumed by a set of variables
         * 
         * @param c the reason clause for a variable in the learnt clause
         * @param inLearnt true if a variable appears in the learnt clause
         * @return true if the reason clause is subsumed, false otherwise
         * 
         * @note this is a helper method for @code{simplifyClauseBasic}
         */
        bool reasonSubsumed(const Clause& c, vec<char>& inLearnt);

        /**
         * @brief Simplify a learnt clause
         * 
         * @param learntClause the learnt clause to modify.
         */
        void simplifyClauseBasic(vec<Lit>& learntClause);

        /**
         * @brief Learn a clause according to the first UIP learning scheme.
         * @param confl the conflict clause
         * @param learntClause the output learnt clause. Assumed to be empty initially.
         */
        void getFirstUIPClause(CRef confl, vec<Lit>& learntClause);

        /**
         * @brief Simplify a learnt clause
         * 
         * @param learntClause the learnt clause to modify.
         */
        void simplifyClause(vec<Lit>& learntClause);
        
        /**
         * @brief Enforce the watcher invariant for a learnt clause. Assumes that the variable with
         * the highest decision level is at index 0. Moves the variable with the next-highest
         * decision level to index 1.
         * 
         * @param learntClause the learnt clause to modify.
         */
        void enforceWatcherInvariant(vec<Lit>& learntClause);

    public:
        //////////////////
        // CONSTRUCTORS //
        //////////////////

        ConflictAnalyzer(Solver* s);
        ~ConflictAnalyzer() = default;

        ////////////////
        // PUBLIC API //
        ////////////////

        void newVar(Var v);

        void analyze     (CRef confl, vec<Lit>& out_learnt, int& out_btlevel); // (bt = backtrack)
        void analyzeFinal(Lit p, vec<Lit>& out_conflict);                      // COULD THIS BE IMPLEMENTED BY THE ORDINARIY "analyze" BY SOME REASONABLE GENERALIZATION?
    };

    inline void ConflictAnalyzer::newVar(Var v) {
    }
}

#endif