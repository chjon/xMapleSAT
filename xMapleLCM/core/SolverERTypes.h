/*******************************************************************************************[SolverERTypes.h]
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

#ifndef Minisat_SolverERTypes_h
#define Minisat_SolverERTypes_h

#include <functional>
#include <tr1/unordered_set>
#include <utility>
#include <vector>
#include <core/SolverTypes.h>

#ifndef ER_ENABLE_GLUCOSER
    #define ER_ENABLE_GLUCOSER false
#endif
#ifndef ER_USER_ADD_SUBEXPR_SET_INTERSECTION
    #define ER_USER_ADD_SUBEXPR_SET_INTERSECTION true
#endif
#ifndef EXTENSION_SUBSTITUTION
    #define EXTENSION_SUBSTITUTION true
#endif
#ifndef EXTENSION_FORCE_BRANCHING
    #define EXTENSION_FORCE_BRANCHING false
#endif
#ifndef ER_PRIORITIZE_EXTVAR
    #define ER_PRIORITIZE_EXTVAR true
#endif

// Define heuristic for filtering clauses before clause selection
#define ER_FILTER_HEURISTIC_NONE     0 // Consider all clauses
#define ER_FILTER_HEURISTIC_RANGE    1 // Consider clauses whose widths are in a certain range 
#define ER_FILTER_HEURISTIC_LONGEST  2 // Consider the longest clauses
#define ER_FILTER_HEURISTIC_LBD      3 // Consider clauses whose LBDs are in a certain range
#define ER_FILTER_HEURISTIC_GLUCOSER 4 // Consider the most recently learnt clauses
#ifndef ER_USER_FILTER_HEURISTIC
    #define ER_USER_FILTER_HEURISTIC ER_FILTER_HEURISTIC_GLUCOSER
#endif

// Define heuristic for selecting clauses
#define ER_SELECT_HEURISTIC_ALL       0 // Consider all clauses
#define ER_SELECT_HEURISTIC_ACTIVITY  1 // Select most active clauses
#define ER_SELECT_HEURISTIC_ACTIVITY2 2 // Select most active clauses using quickselect
#define ER_SELECT_HEURISTIC_GLUCOSER  3 // Only consider the previous two learnt clauses
#ifndef ER_USER_SELECT_HEURISTIC
    #define ER_USER_SELECT_HEURISTIC ER_SELECT_HEURISTIC_GLUCOSER
#endif

// Define heuristic for replacing extension definitions in clauses
#define ER_SUBSTITUTE_HEURISTIC_NONE  0x0 // Don't substitute
#define ER_SUBSTITUTE_HEURISTIC_WIDTH 0x1 // Consider clauses within a clause width range
#define ER_SUBSTITUTE_HEURISTIC_LBD   0x2 // Consider clauses within an LBD range
#ifndef ER_USER_SUBSTITUTE_HEURISTIC
    #define ER_USER_SUBSTITUTE_HEURISTIC (ER_SUBSTITUTE_HEURISTIC_WIDTH | ER_SUBSTITUTE_HEURISTIC_LBD)
#endif

// Define heuristics for adding extension definitions
#define ER_ADD_HEURISTIC_NONE     0 // Do not add extension variables
#define ER_ADD_HEURISTIC_RANDOM   1 // Add extension variables by selecting random pairs of literals
#define ER_ADD_HEURISTIC_SUBEXPR  2 // Add extension variables by selecting the most common pairs of literals
#define ER_ADD_HEURISTIC_GLUCOSER 3 // Add extension variables according to the scheme prescribed by GlucosER
#ifndef ER_USER_ADD_HEURISTIC
    #define ER_USER_ADD_HEURISTIC ER_ADD_HEURISTIC_GLUCOSER
#endif

// Define heuristics for deleting extension variables
#define ER_DELETE_HEURISTIC_NONE      0 // Do not delete extension variables
#define ER_DELETE_HEURISTIC_ALL       1 // Delete all extension variables
#define ER_DELETE_HEURISTIC_ACTIVITY  2 // Delete low-activity extension variables based on a constant activity threshold
#define ER_DELETE_HEURISTIC_ACTIVITY2 3 // Delete low-activity extension variables based on a proportional activity threshold
#ifndef ER_USER_DELETE_HEURISTIC
    #define ER_USER_DELETE_HEURISTIC ER_DELETE_HEURISTIC_NONE
#endif

// Define heuristics for location to generate variables
#define ER_GEN_LOCATION_NONE           0 // Do not generate extension variables
#define ER_GEN_LOCATION_AFTER_RESTART  1 // Generate extension variables after a restart
#define ER_GEN_LOCATION_AFTER_CONFLICT 2 // Generate extension variables after a conflict
#ifndef ER_USER_GEN_LOCATION
    #define ER_USER_GEN_LOCATION ER_GEN_LOCATION_AFTER_RESTART
#endif

// Define heuristics for location to add variables
#define ER_ADD_LOCATION_NONE           0 // Do not add extension variables
#define ER_ADD_LOCATION_AFTER_RESTART  1 // Add extension variables after a restart
#define ER_ADD_LOCATION_AFTER_CONFLICT 2 // Add extension variables after a conflict
#ifndef ER_USER_ADD_LOCATION
    #define ER_USER_ADD_LOCATION ER_GEN_LOCATION_AFTER_RESTART
#endif

namespace Minisat {

using LitPair = std::pair<Lit, Lit>;
using LitSet  = std::tr1::unordered_set<Lit>;

inline LitPair mkLitPair(Lit a, Lit b) {
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
    std::vector< std::vector<Lit> > additionalClauses;

    friend std::ostream &operator<<(std::ostream &os, const ExtDef& def) { 
        os << "{" << def.x << "=" << def.a << "v" << def.b;
        if (def.additionalClauses.size()) {
            os << ":[";
            for (const std::vector<Lit>& clause : def.additionalClauses) {
                os << "[";
                if (clause.size()) os << clause[0];
                for (Lit l : clause) os << "," << l; 
                os << "]";
            }
            os << "]";
        }
        return os << "}";
    }
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

// Extension Variable Deletion

/**
 * @brief A user-defined function to set up for determining whether to delete a variable
 */
using DeletionPredicateSetup = std::function<void()>;

/**
 * @brief A user-defined predicate for determining whether to delete a variable
 * 
 * @param x the candidate variable to delete
 * @return true if the solver should delete the variable
 * @return false otherwise
 */
using DeletionPredicate = std::function<bool(Var)>;

}

#endif