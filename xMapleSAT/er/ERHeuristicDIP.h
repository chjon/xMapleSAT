/********************************************************************************[ERHeuristicDIP.h]
xMaple* -- Copyright (c) 2022, Jonathan Chung, Vijay Ganesh, Sam Buss

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

#ifndef Minisat_ERManager_h
#define Minisat_ERManager_h

#include <vector>
#include "core/AssignmentTrail.h"
#include "core/SolverTypes.h"

namespace Minisat {
    class ERHeuristicDIP {
    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // MEMBER VARIABLES

        /// @brief Whether variables have been seen while exploring the conflict graph
        vec<bool> seen;

        /// @brief Literals that participate in the conflict graph, in reverse topological order
        vec<Lit> conflictLits;

        /// @brief Remapping variables in the conflict graph to the range [0, N]
        vec<int> remappedVariables;

        /// @brief The predecessors of a variable in the conflict graph, indexed according to
        /// @code{predecessorIndex}
        vec<int> predecessors;

        /**
         * @brief Indices into @code{predecessors}.
         * @details For i < N-1, if p = predecessorIndex[i], then predecessors[p+k] is the k-th
         * predecessor of i. predecessorIndex[N-1] is the number of entries in predecessors.
         * 
         */
        vec<int> predecessorIndex;

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // SOLVER REFERENCES
        AssignmentTrail& assignmentTrail;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CONSTRUCTORS / DESTRUCTORS

        ERHeuristicDIP(Solver& s);
        ~ERHeuristicDIP() = default;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // VARIABLE DEFINITION HEURISTIC API

        /**
         * @brief A user-defined method for generating extension variable definitions given a list
         * of selected clauses
         * 
         * @param extVarDefBuffer the buffer in which to output extension variable definitions
         * @param maxNumNewVars the preferred number of extension variables to define.
         */
        void generateDefinitions(std::vector<ExtDef>& extVarDefBuffer, unsigned int maxNumNewVars);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // EVENT HANDLERS

        /**
         * @brief Set up internal data structures for a new variable
         * 
         * @param v the variable to register
         */
        void newVar(Var v);

        /**
         * @brief Update internal data structures when the solver encounters a conflict.
         * @details Clears data structures to prepare to populate @code{predecessor} and
         * @code{predecessorIndex}.
         * 
         */
        void handleEventConflicted(void);

        /**
         * @brief Update internal data structures when a literal appears in the conflict graph.
         * @details Adds literals to a stack to prepare to populate @code{predecessor} and
         * @code{predecessorIndex}.
         * 
         * @param l The literal that appears in the conflict graph
         * 
         * @note Assumes that literals appear in reverse topological order with respect to the
         * conflict graph.
         */
        void handleEventLitInConflictGraph(Lit l);

        /**
         * @brief Update internal data structures when a clause is learned.
         * @details Populates the @code{predecessor} and @code{predecessorIndex} arrays.
         * 
         * @param c The learnt clause
         * 
         * @pre @code{conflictLits} contains all literals in the conflict graph between the
         * conflict and the first UIP literal in reverse topological order.
         */
        void handleEventLearntClause(const Clause& c);
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE METHODS

    inline void ERHeuristicDIP::newVar(Var v) {
        seen.push(false);
        remap.push(0);
    }

    inline void ERHeuristicDIP::handleEventLitInConflictGraph(Lit l) {
        conflictLits.push(l);
    }
}

#endif