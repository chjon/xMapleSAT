/*******************************************************************************[ERHeuristicDIP.cc]
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

#include "er/ERHeuristicDIP.h"
#include "er/ERTypes.h"
#include "er/TwoVertexBottlenecks.h"
#include "core/Solver.h"

using namespace Minisat;

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

ERHeuristicDIP::ERHeuristicDIP(Solver& s)
    : assignmentTrail(s.assignmentTrail)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// VARIABLE DEFINITION HEURISTIC API

void ERHeuristicDIP::generateDefinitions(
    std::vector<ExtDef>& extVarDefBuffer,
    unsigned int maxNumNewVars // Ignoring this parameter for now
) {
    // Compute DIP using predecessors
    TwoVertexBottlenecks dipCollector(
        user_extDefHeuristic_dip_stack.size(),
        user_extDefHeuristic_dip_predecessors,
        user_extDefHeuristic_dip_predecessorIndex
    );
    dipCollector.computeDIPs();

    // Exit if there is no DIP 
    if (dipCollector.NumPairGroups() == 0) return;
    
    // Get DIP literals
    const Var z = assignmentTrail.nVars() + extVarDefBuffer.size();
    const int dipIndexA = dipCollector.GroupMemberLastLeft();
    const int dipIndexB = dipCollector.GroupMemberLastRight();
    const Lit a = getLitFromIndex(stack, dipIndexA);
    const Lit b = getLitFromIndex(stack, dipIndexB);

    // Add definition
    // Note: we want (z' = a AND b), but we define vars using disjunctions. Then we can define
    // (z = ~z'), so (z = ~a OR ~b)
    extVarDefBuffer.push(ExtDef{ mkLit(z), ~a, ~b, std::vector< std::vector<int> >() });

    // Add helper clauses
    std::vector<Lit> helperClause1 = getLowerLevelLitsBeforeDIP(dipIndexA, dipIndexB); // 'D'
    std::vector<Lit> helperClause2 = getLowerLevelLitsAfterDIP(dipIndexA, dipIndexB);  // 'C'

    const Lit uip = stack.last();
    ExtDef& extDef = extVarDefBuffer[extVarDefBuffer.size() - 1];
    helperClause1.push(~uip, ~z);
    helperClause2.push(z);

    extDef.additionalClauses.push(helperClause1);
    extDef.additionalClauses.push(helperClause2);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// EVENT HANDLERS

void ERHeuristicDIP::handleEventLearntClause(const Clause& c) {
    // Read the stack in reverse to get literals in topological order
    for (int i = stack.size() - 1, numVarsSeen = 0; i >= 0; i--, numVarsSeen++) {
        const Var v = var(stack[i]);
        
        // Remap conflict graph variables to consecutive indices starting from 0
        // TODO: check whether the DIP algo wants predecessors to contain vars at earlier levels
        remappedVariables[v] = numVarsSeen;
        predecessorIndex[remappedVariables[v]] = predecessors.size();

        // Don't store predecessors for the first UIP or for variables at earlier levels
        if (numVarsSeen == 0 ||
            assignmentTrail.level(v) < assignmentTrail.currentDecisionLevel()
        ) {
            continue;
        }

        // Store predecessors (remapped)
        const Clause& reasonClause = ca[assignmentTrail.reason(v)];
        for (int i = 0; i < reasonClause.size(); i++) {
            const Var x = var(reasonClause[i]);
            predecessors.push(remappedVariables[x]);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

/**
 * @brief Get the literal in the conflict graph from the index returned by the DIP algorithm
 * 
 * @param stack the literals in the conflict graph in reverse topological order
 * @param index the index for which to get the literal
 * @return the literal corresponding to the requested index
 */
static inline Lit getLitFromIndex(const vec<Lit>& stack, int index) {
    return stack[stack.size() - index - 1];
}

/**
 * @brief Get the literals from earlier levels that participate in the conflict graph before the DIP
 * 
 * @param stack the literals in the conflict graph in reverse topological order
 * @param indexA the index of the first DIP literal
 * @param indexB the index of the second DIP literal
 * @return a list of literals from earlier levels that participate in the conflict graph before the DIP
 */
static inline std::vector<Lit> getLowerLevelLitsBeforeDIP(const vec<Lit>& stack, int indexA, int indexB);

/**
 * @brief Get the literals from earlier levels that participate in the conflict graph after the DIP
 * 
 * @param stack the literals in the conflict graph in reverse topological order
 * @param indexA the index of the first DIP literal
 * @param indexB the index of the second DIP literal
 * @return a list of literals from earlier levels that participate in the conflict graph after the DIP
 */
static inline std::vector<Lit> getLowerLevelLitsAfterDIP(const vec<Lit>& stack, int indexA, int indexB);