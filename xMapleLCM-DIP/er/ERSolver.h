/**************************************************************************************[ERSolver.h]
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

#ifndef Minisat_ERSolver_h
#define Minisat_ERSolver_h

#include "core/Solver.h"
#include "er/ERTypes.h"
#include "er/ERManager.h"

namespace Minisat {
    /**
     * @brief This class augments the base solver with extended resolution capabilities.
     * 
     */
    class ERSolver : public Solver {

      // In DIP-learning we might want to learn the clause UIP ^ previous_DL_literals -> DIP
      vec<Lit> learnt_clause_UIP_to_DIP; 

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // SOLVER COMPONENTS

        /// @brief The extended resolution component
        ERManager erManager;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CONSTRUCTORS

        /**
         * @brief Construct a new ERSolver object
         * 
         */
        ERSolver();

        /**
         * @brief Destroy the ERSolver object
         * 
         */
        virtual ~ERSolver() = default;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATE MODIFICATION

        /**
         * @brief Relocate all clauses
         * 
         * @param to the ClauseAllocator to relocate to
         */
        virtual void relocAll(ClauseAllocator& to);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // PROBLEM SPECIFICATION

        /**
         * @brief Creates a new SAT variable in the solver
         * 
         * @param polarity The preferred initial polarity of the new variable
         * @param dvar true iff the new variable should be used as a decision variable (NOTE! This has
         * effects on the meaning of a SATISFIABLE result).
         * @return the ID of the new variable
         */
        virtual Var newVar (bool polarity = true, bool dvar = true);

    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CORE SEARCH FUNCTIONS

        /**
         * @brief Search for a given number of conflicts
         * 
         * @param nof_conflicts The number of conflicts to search for before exiting the function
         * @return 'l_True' if a partial assigment that is consistent with respect to the clause set
         * is found. If all variables are decision variables, this means that the clause set is
         * satisfiable. 'l_False' if the clause set is unsatisfiable. 'l_Undef' if the bound on number
         * of conflicts is reached.
         */
        virtual lbool search(int& nof_conflicts);
        
        /**
         * @brief Main solve method (assumptions given in 'assumptions')
         * 
         * @return l_True if the formula is satisfiable, l_False if the formula is unsatisfiable, or
         * l_Undef if the satisfiability of the formula is unknown.
         */
        virtual lbool solve_(void);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // EVENT HANDLERS

        /**
         * @brief Update data structures when a clause is deleted
         * 
         * @param c the clause that was deleted
         * @param cr the reference to the clause that was deleted
         */
        virtual void handleEventClauseDeleted(const Clause& c, CRef cr);
      
    };

  
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE FUNCTIONS

    ///////////////////////
    // STATE MODIFICATION

    inline void ERSolver::relocAll(ClauseAllocator& to) {
        Solver::relocAll(to);
        erManager.relocAll(to);
    }

    //////////////////////////
    // PROBLEM SPECIFICATION

    inline Var ERSolver::newVar(bool sign, bool dvar) {
        Var v = Solver::newVar(sign, dvar);
        erManager.newVar(v);
        return v;
    }

    ///////////////////
    // EVENT HANDLERS

    inline void ERSolver::handleEventClauseDeleted(const Clause& c, CRef cr) {
        Solver::handleEventClauseDeleted(c, cr);
        erManager.handleEventClauseDeleted(cr);
    }
}

#endif
