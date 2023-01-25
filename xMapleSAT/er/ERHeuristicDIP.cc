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

#include <vector>
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
    // Set up data structures for TwoVertexBottlenecks
    dipComputationSetup(conflictLits, predecessors, predecessorIndex, seen, assignmentTrail, ca, confl);

    // Compute DIP using predecessors
    TwoVertexBottlenecks dipCollector(
        predecessorIndex.size(),
        (int*)predecessors,
        (int*)predecessorIndex
    );
    dipCollector.computeDIPs();

    // Exit if there is no DIP 
    if (dipCollector.NumPairGroups() == 0) return;
    
    // Get DIP literals
    const Var z = assignmentTrail.nVars() + extVarDefBuffer.size();
    const int dipIndexA = dipCollector.GroupMemberLastLeft();
    const int dipIndexB = dipCollector.GroupMemberLastRight();
    const Lit a = getLitFromIndex(conflictLits, dipIndexA);
    const Lit b = getLitFromIndex(conflictLits, dipIndexB);

    // Add definition
    // Note: we want (z' = a AND b), but we define vars using disjunctions. Then we can define
    // (z = ~z'), so (z = ~a OR ~b)
    extVarDefBuffer.push(ExtDef{ mkLit(z), ~a, ~b, std::vector< std::vector<int> >() });

    // Add helper clauses
    std::vector<Lit> helperClause1 = getLowerLevelLitsAfterDIP(litsAfterDIP, dipA, dipB, confl, assignmentTrail, ca); // 'D'
    std::vector<Lit> helperClause2 = getLowerLevelLitsBeforeDIP(litsBeforeDIP, dipA, dipB, assignmentTrail, ca);      // 'C'

    const Lit uip = conflictLits.last();
    ExtDef& extDef = extVarDefBuffer[extVarDefBuffer.size() - 1];
    helperClause1.push(~uip, ~z);
    helperClause2.push(z);

    extDef.additionalClauses.push(helperClause1);
    extDef.additionalClauses.push(helperClause2);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS

/**
 * @brief Get the literals that participate in the conflict graph (up to the first UIP) in reverse
 * topological order.
 * 
 * @details Returns literals at the current decision level and literals at lower decision levels
 * that participate in propagating literals at the current decision level 
 * 
 * @param conflictLits 
 * @param seen 
 * @param assignmentTrail 
 * @param ca 
 * @param confl 
 */
static void getFirstUIPConflictLits(
    vec<Lit>& conflictLits,
    vec<bool>& seen,
    AssignmentTrail& assignmentTrail,
    ClauseAllocator& ca,
    CRef confl
) {
    // Initialize local data structures
    Lit p = lit_Undef;
    int pathC = 0;
    int index = assignmentTrail.nAssigns() - 1;
    
    // Iterate through every clause that participates in the conflict graph, backtracking until the first UIP
    do {
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

        // Iterate through every literal that participates in the conflict graph
        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];
            if (seen[var(q)] || assignmentTrail.level(var(q)) == 0) continue;

            // Mark variable as seen
            seen[var(q)] = true;

            // Increment the number of paths if the variable is assigned at the decision level that resulted in conflict
            if (assignmentTrail.level(var(q)) == assignmentTrail.decisionLevel()) {
                pathC++;

                // Add variable to list of conflict literals
                conflictLits.push(q);
            }
        }
        
        // Select next clause to look at:
        while (!seen[var(assignmentTrail[index--])]);
        p     = assignmentTrail[index+1];
        confl = assignmentTrail.reason(var(p));

        // Mark variable as unseen: it is either at or after the first UIP
        pathC--;
    } while (pathC > 0);
}

/**
 * @brief Map conflict variables to the range [0,N-1] in topological order
 * 
 * @param conflictLits 
 * @param predecessors 
 * @param predecessorIndex 
 * @param assignmentTrail 
 * @param ca 
 */
static void remapConflictLits(
    vec<Lit>& conflictLits,
    vec<int>& predecessors,
    vec<int>& predecessorIndex,
    vec<int>& remappedVariables,
    AssignmentTrail& assignmentTrail,
    ClauseAllocator& ca
) {
    for (int i = conflictLits.size() - 1; i >= 0; i--) {
        const Var v = var(conflictLits[i]);

        // Store index and remap current variable
        predecessorIndex.push(predecessors.size());
        remappedVariables[v] = predecessorIndex.size();

        // Store remapped predecessors
        CRef reason = assignmentTrail.reason(v);
        if (i = conflictLits.size() - 1 || reason == CRef_Undef) continue;
        const Clause& reasonClause = ca[reason];
        for (int j = 0; j < reasonClause.size(); j++) {
            if (assignmentTrail.level() != assignmentTrail.decisionLevel()) continue;
            const Var x = var(reasonClause[j]);
            predecessors.push(remappedVariables[x]);
        }
    }
}

static void dipComputationSetup(
    vec<Lit>& conflictLits,
    vec<int>& predecessors,
    vec<int>& predecessorIndex,
    vec<int>& remappedVariables,
    vec<bool>& seen,
    AssignmentTrail& assignmentTrail,
    ClauseAllocator& ca,
    CRef confl
) {
    // Clear data structures
    conflictLits.clear();
    predecessors.clear();
    predecessorIndex.clear();

    // Get conflict literals in reverse topological order
    getFirstUIPConflictLits(conflictLits, seen, assignmentTrail, ca, confl);

    // Get conflict variables in topological order and remap them to the range [0,N-1]
    remapConflictLits(conflictLits, predecessors, predecessorIndex, remappedVariables, assignmentTrail, ca);

    // Clean up
    for (int i = 0; i < conflictLits.size(); i++)
        seen[var(conflictLits[i])] = false;
}

/**
 * @brief Get the literal in the conflict graph from the index returned by the DIP algorithm
 * 
 * @param conflictLits the literals in the conflict graph in reverse topological order
 * @param index the index for which to get the literal
 * @return the literal corresponding to the requested index
 */
static inline Lit getLitFromIndex(const vec<Lit>& conflictLits, int index) {
    return conflictLits[conflictLits.size() - index - 1];
}

/**
 * @brief Get the literals from earlier levels that participate in the conflict graph after the DIP
 * 
 * @param dipA the first DIP literal
 * @param dipB the second DIP literal
 * @return a list of literals from earlier levels that participate in the conflict graph after the DIP
 */
static inline std::vector<Lit> getLowerLevelLitsAfterDIP(
    Var dipA,
    Var dipB,
    CRef confl,
    AssignmentTrail& assignmentTrail,
    ClauseAllocator& ca
) {
    // Initialize local data structures
    Lit p = lit_Undef;
    int pathC = 0;
    int dipSeen = 0;
    int index = 0;
    vec<Var> bfsQueue;

    std::vector<Lit> litsAfterDIP;

    // Iterate through variables that participate in the conflict graph, backtracking until the DIP
    do {
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

        // Iterate through every literal that participates in the conflict graph
        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];
            if (seen[var(q)] || assignmentTrail.level(var(q)) == 0) continue;

            // Mark variable as seen
            seen[var(q)] = true;
            bfsQueue.push(var(q));

            // Add variable to litsAfterDIP
            if (assignmentTrail.level(var(q)) < assignmentTrail.decisionLevel())
                litsAfterDIP.push_back(q);

            // Increment the number of paths if the variable is assigned at the decision level that resulted in conflict
            if (assignmentTrail.level(var(q)) == assignmentTrail.decisionLevel())
                pathC++;
        }
        
        // Select next clause to look at:
        do {
            p = bfsQueue[index++];
            confl = assignmentTrail.reason(var(p));
            seen[var(p)] = false;
            pathC--;

            if (var(p) == dip1 || var(p) == dip2) {
                dipSeen++;
            } else {
                break;
            }
        } while (pathC > 0);
    } while (dipSeen < 2 && pathC > 0); // Shouldn't need to check pathC if the DIP computation algorithm worked 

    assert(pathC == 2);

    // Clean up
    for (int i = 0; i < litsAfterDIP.size(); i++)
        seen[litsAfterDIP[i]] = false;

    return litsAfterDIP;
}

/**
 * @brief Get the literals from earlier levels that participate in the conflict graph before the DIP
 * 
 * @param dipA the first DIP literal
 * @param dipB the second DIP literal
 * @return a list of literals from earlier levels that participate in the conflict graph before the DIP
 */
static inline std::vector<Lit> getLowerLevelLitsBeforeDIP(
    Var dipA,
    Var dipB,
    AssignmentTrail& assignmentTrail,
    ClauseAllocator& ca
) {
    // Set up data structures to iterate backward from DIP
    CRef confl = reason(dip1);
    Lit p = mkLit(dip1);
    int pathC = 0;
    int index = 0;
    vec<Var> bfsQueue;
    bfsQueue.push(dip2);

    std::vector<Lit> litsBeforeDIP;

    // Iterate through variables that participate in the conflict graph, backtracking until the UIP
    do {
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

        // Iterate through every literal that participates in the conflict graph
        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];
            if (seen[var(q)] || assignmentTrail.level(var(q)) == 0) continue;

            // Mark variable as seen
            seen[var(q)] = true;
            bfsQueue.push(var(q));

            // Add variable to litsBeforeDIP
            if (assignmentTrail.level(var(q)) < assignmentTrail.decisionLevel())
                litsBeforeDIP.push_back(q);

            // Increment the number of paths if the variable is assigned at the decision level that resulted in conflict
            if (assignmentTrail.level(var(q)) == assignmentTrail.decisionLevel())
                pathC++;
        }
        
        // Select next clause to look at:
        p = bfsQueue[index++];
        seen[var(p)] = false;
        confl = assignmentTrail.reason(var(p));
        pathC--;
    } while (pathC > 0);

    // Clean up
    for (int i = 0; i < litsBeforeDIP.size(); i++)
        seen[litsBeforeDIP[i]] = false;

    return litsBeforeDIP;
}