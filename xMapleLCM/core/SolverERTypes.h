#ifndef Minisat_SolverERTypes_h
#define Minisat_SolverERTypes_h

#include <functional>
#include <utility>
#include <vector>
#include <core/SolverTypes.h>

namespace Minisat {

inline std::pair<Lit, Lit> mkLitPair(Lit a, Lit b) {
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

// Clause Selection

/**
 * @brief A user-defined method for filtering clauses for clause selection 
 * 
 * @param candidate the clause to check
 * 
 * @return true if the clause satisfies the predicate
 * @return false otherwise
 */
using FilterPredicate = std::function<bool(CRef)>;

/**
 * @brief A user-defined method for selecting clauses for defining extension variables
 * 
 * @param output the output list of selected clauses
 * @param input the input list from which to select clauses
 * @param numClauses the number of clauses to select
 */
using SelectionHeuristic = std::function<void(std::vector<CRef>&, const std::vector<CRef>&, unsigned int)>;

// Extension Variable Definition

/**
 * @brief A structure representing an extension variable definition
 */
struct ExtDef {
    // The positive literal representing the extension variable
    Lit x;
    
    // The two basis literals such that (x <=> a v b) 
    Lit a, b;
    
    // Clauses to learn in addition to the three defining clauses encoding (x <=> a v b) 
    std::vector<std::vector<Lit>> additionalClauses;
};

/**
 * @brief A user-defined method for generating extension variable definitions given a list of selected clauses
 * 
 * @param extVarDefBuffer the buffer in which to output extension variable definitions
 * @param selectedClauses the list of selected clauses from which to generate extension variable definitions
 * @param maxNumNewVars the preferred number of extension variables to define
 */
using ExtDefHeuristic = std::function<void(std::vector<ExtDef>&, const std::vector<CRef>&, unsigned int)>;

// Extension Variable Substitution

/**
 * @brief A user-defined predicate for determining whether to attempt substitution into a clause
 * 
 * @param clause the vector of literals to check 
 * @return true if the solver should try substituting into the clause
 * @return false otherwise
 */
using SubstitutionPredicate = std::function<bool(vec<Lit>&)>;

}

#endif