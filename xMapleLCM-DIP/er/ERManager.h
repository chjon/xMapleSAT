/*************************************************************************************[ERManager.h]
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

#include <functional>
#include <initializer_list>
#include <map>
#include <tuple>
#include <utility>
#include <vector>
#include <sys/resource.h>
#include <tr1/unordered_set>
#include <tr1/unordered_map>

#include "core/AssignmentTrail.h"
#include "core/SolverTypes.h"
#include "er/ERTypes.h"
#include "mtl/Vec.h"
#include "mtl/ExtDefMap.h"

using namespace std;

// Making some internal methods visible for testing
#ifdef TESTING
#define protected public
#endif

namespace Minisat {

// Forward declarations
class BranchingHeuristicManager;
class ClauseDatabase;
class PropagationQueue;
class RandomNumberGenerator;
class UnitPropagator;

class ERManager {
public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // HELPER TYPES

    /**
     * @brief The overall definition heuristic
     * @note This is used to choose separate data structures for different heuristics
     */
    enum HeuristicType : int {
        DEFAULT = 0,
        LER = 1
    };

private:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // SOLVER REFERENCES

    AssignmentTrail& assignmentTrail;
    BranchingHeuristicManager& branchingHeuristicManager;
    ClauseDatabase& clauseDatabase;
    PropagationQueue& propagationQueue;
    RandomNumberGenerator& randomNumberGenerator;
    UnitPropagator& unitPropagator;
    ClauseAllocator& ca;
    Solver& solver;

private:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // MEMBER VARIABLES

    /// @brief The number of variables in the original formula.
    /// @note This value is used to quickly check whether a variable is an extension variable
    int originalNumVars;

    /// @brief Map from extension variables to a list of their extension definition clauses.
    std::tr1::unordered_map<Var, std::vector<CRef> > extDefs;

    /// @brief Map from variables to their extension level
    vec<unsigned int> extensionLevel;

    /// @brief Map from extension variables to the pair of literals they represent
    /// @note This is used for extension variable substitution
    ExtDefMap<Lit> xdm;

    /// @brief Stores the last time extension variables were added.
    /// @note This is used to determine when to introduce new extension variables
    uint64_t prevExtensionConflict;

    /// @brief Stores the last time extension variables were deleted.
    /// @note This is used to determine when to delete extension variables
    uint64_t prevDelExtVarConflict;
    
    std::vector<CRef> m_filteredClauses;
    std::vector<CRef> m_selectedClauses;
    std::vector<ExtDef> m_extVarDefBuffer;

    /// @brief Keep track of clauses which have been removed so that they can be removed from the
    /// buffers above.
    std::tr1::unordered_set<CRef> m_deletedClauses;

private:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // DATA STRUCTURES FOR USER-DEFINED HEURISTICS
    
    // Local Extended Resolution (LER)

    std::vector<CRef> m_filteredClauses_ler;
    std::vector<CRef> m_selectedClauses_ler;
    std::vector<ExtDef> m_extVarDefBuffer_ler;

private:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TEMPORARY DATA STRUCTURES

    // These are declared here to avoid repeated memory allocation.
    // These data structures should only be used locally.
    
    vec<Lit> tmp_vec;
    LitSet tmp_set;

private:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // HEURISTIC FUNCTIONS

    FilterPredicate        user_extFilPredicate;
    SelectionHeuristic     user_extSelHeuristic;
    ExtDefHeuristic        user_extDefHeuristic;
    SubstitutionPredicate  user_extSubPredicate;
    DeletionPredicateSetup user_extDelPredicateSetup;
    DeletionPredicate      user_extDelPredicate;

protected:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // PARAMETERS

    /// @brief Number of conflicts to wait before trying to introduce an extension variable.
    /// (default 2000)
    int ext_freq;

    /// @brief Number of clauses to consider when introducing extension variables. (default 100)
    int ext_window;

    /// @brief Maximum number of extension variables to introduce at once. (default 1)
    int ext_max_intro;

    /// @brief The fraction of maximum activity that should be given to new variables (default 1.0)
    double ext_prio_act;

    /// @brief Preferred sign for new variables (default true (negated))
    bool ext_pref_sign;

    /// @brief Minimum clause width to consider when selecting clauses
    int ext_min_width;

    /// @brief Maximum clause width to consider when selecting clauses
    int ext_max_width;

    /// @brief Maximum number of clauses after the filter step
    int ext_filter_num;

    /// @brief Minimum width of clauses to substitute into
    int ext_sub_min_width;

    /// @brief Maximum width of clauses to substitute into
    int ext_sub_max_width;

    /// @brief Minimum LBD of clauses to substitute into
    int ext_min_lbd;

    /// @brief Maximum LBD of clauses to substitute into
    int ext_max_lbd;

    /// @brief Activity threshold for deleting clauses
    double ext_act_threshold;

    /// @brief Number of conflicts to wait before trying to delete extension variables
    int ext_del_freq;

    /// @brief Threshold activity for variable deletion
    double m_threshold_activity;

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // STATISTICS

    /// @brief The total number of extension variables introduced by the solver
    mutable uint64_t total_ext_vars;

    /// @brief The total number of extension variables which the solver attempted to be deleted
    mutable uint64_t tried_del_ext_vars;

    /// @brief The total number of extension variables which the solver successfully deleted
    mutable uint64_t deleted_ext_vars;

    /// @brief The maximum number of extension variables present in the solver at any one time
    mutable uint64_t max_ext_vars;

    /// @brief The total number of clauses containing extension variables which participated in a
    /// conflict
    mutable uint64_t conflict_extclauses;
    
    /// @brief The total number of learnt clauses containing extension variables 
    mutable uint64_t learnt_extclauses;

    /// @brief The total number of times the solver branched on an extension variable
    mutable uint64_t branchOnExt;

    ///////////////////////////////////////////
    // Measuring extended resolution overhead

    mutable struct rusage ext_timer_start, ext_timer_end;
    mutable struct rusage ext_sel_overhead; // Overhead for selecting clauses for adding extension variables
    mutable struct rusage ext_add_overhead; // Overhead for adding extension variables
    mutable struct rusage ext_delC_overhead; // Overhead for deleting clauses containing extension variables
    mutable struct rusage ext_delV_overhead; // Overhead for deleting extension variables
    mutable struct rusage ext_sub_overhead; // Overhead for substituting disjunctions containing extension variables
    mutable struct rusage ext_stat_overhead; // Overhead for measuring statistics

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // CONSTRUCTORS

    /**
     * @brief Construct a new ERManager object
     * 
     * @param s a reference to the main solver object
     */
    ERManager(Solver& s);

    /**
     * @brief Destroy the ERManager object
     * 
     */
    ~ERManager() = default;

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // ACCESSORS

    /**
     * @brief Determine whether a variable is an extension variable
     * 
     * @param x the variable to check
     * @return true if the variable is an extension variable, regardless of whether it has been
     * deleted
     * @return false otherwise
     */
    inline bool isExtVar(Var x) const;

    /**
     * @brief Determine whether a variable is an extension variable and has not been deleted
     * 
     * @param x the variable to check
     * @return true if the variable is an extension variable and has not been deleted
     * @return false otherwise
     */
    inline bool isCurrentExtVar(Var x) const;

    /**
     * @brief Determine whether a pair of literals can be used as the basis literals for a new
     * extension variable
     * 
     * @param a the first literal
     * @param b the second literal
     * @param generatedPairs the set of literal pairs queued to be added to the solver
     * @return false if the literals refer to the same variable, if the literals have been set at
     * the root level, or the pair of literals has already been added or queued to be added to the
     * solver
     * @return true otherwise
     * 
     * @note This method is intended to be used by variable definition generation heuristics
     */
    bool isValidDefPair(
        Lit a,
        Lit b,
        const std::tr1::unordered_set<LitPair>& generatedPairs
    ) const;

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // STATE MODIFICATION

    /**
     * @brief Set up data structures when the solver starts
     * 
     */
    inline void init(void);

    /**
     * @brief Update data structures to allocate enough memory when a new variable is added
     * 
     * @param v the variable to register
     */
    inline void newVar(Var v);

    /**
     * @brief Relocate CRefs to new ClauseAllocator
     * 
     * @param to The ClauseAllocator into which to reloc 
     */
    void relocAll(ClauseAllocator& to);

    /**
     * @brief Remove extension definition clauses that have already been satisfied
     */
    void removeSatisfied();

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // CLAUSE SELECTION

    /**
     * @brief Filter clauses for clause selection
     * 
     * @param candidate the clause to check
     * @param filterPredicate a method for selecting clauses for defining extension variables
     * @param heuristicType the overall heuristic type
     * 
     * @note Saves outputs in @code{m_filteredClauses}
     */
    void filterIncremental(
        const CRef candidate,
        FilterPredicate& filterPredicate,
        HeuristicType heuristicType = HeuristicType::DEFAULT
    );

    /**
     * @brief Filter clauses for clause selection
     * 
     * @param candidate the clause to check
     * @param heuristicType the overall heuristic type
     * 
     * @note Saves outputs in @code{m_filteredClauses}
     */
    void filterIncremental(const CRef candidate, HeuristicType heuristicType = HeuristicType::DEFAULT);

    /**
     * @brief Select clauses for defining extension variables
     * 
     * @param selectionHeuristic a method for selecting clauses for defining extension variables
     * @param heuristicType the overall heuristic type
     * @param numKeepFiltered the number of filtered clauses to keep between function calls
     * 
     * @note Takes input from @code{m_filteredClauses}
     * @note Saves outputs in @code{m_selectedClauses}
     */
    void selectClauses(
        SelectionHeuristic& selectionHeuristic,
        HeuristicType heuristicType = HeuristicType::DEFAULT,
        unsigned int numKeepFiltered = 0
    );

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // EXTENSION VARIABLE DEFINITION

    /**
     * @brief Generate extension variable definitions
     * 
     * @param extDefHeuristic a method for generating extension variable definitions given a list
     * of selected clauses
     * @param heuristicType the overall heuristic type
     */
    void defineExtVars(
        ExtDefHeuristic& extDefHeuristic,
        HeuristicType heuristicType = HeuristicType::DEFAULT
    );
    
    /**
     * @brief Checks whether to generate definitions and then calls @code{selectClauses} and
     * @code{defineExtVars}.
     * 
     * @param conflicts the current total number of conflicts seen by the solver
     */
    inline void checkGenerateDefinitions(uint64_t conflicts);

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // EXTENSION VARIABLE INTRODUCTION

    /**
     * @brief Introduce extension variables into the solver
     * 
     * @param ext_def_db The map of extension variables to a list of their corresponding extension
     * definition clauses. In most cases this should be equal to @code{extDefs}.
     * @param heuristicType the overall heuristic type
     */
    void introduceExtVars(
        std::tr1::unordered_map<Var, std::vector<CRef> >& ext_def_db,
        HeuristicType heuristicType = HeuristicType::DEFAULT
    );

    /**
     * @brief Introduce extension variables into the solver
     * 
     * @param heuristicType the overall heuristic type
     */
    void introduceExtVars(HeuristicType heuristicType = HeuristicType::DEFAULT);

    /**
     * @brief Adds an extension definition clause to the appropriate clause database
     * 
     * @param db The clause database to which to add the clause
     * @param ext_lit The extension variable corresponding to the given clause
     * @param clause The vector of literals to add as a clause
     * 
     * @pre current decision level must be zero
     */
    void addExtDefClause(std::vector<CRef>& db, Lit ext_lit, vec<Lit>& clause);

    /**
     * @brief Adds an extension definition clause to the appropriate clause database
     * 
     * @param db The clause database to which to add the clause
     * @param ext_lit The extension variable corresponding to the given clause
     * @param clause The vector of literals to add as a clause
     * 
     * @pre current decision level must be zero
     */
    void addExtDefClause(
        std::vector<CRef>& db,
        Lit ext_lit,
        const std::initializer_list<Lit>& clause
    );

    /**
     * @brief Adds an extension definition clause to the appropriate clause database
     * 
     * @param db The clause database to which to add the clause
     * @param ext_lit The extension variable corresponding to the given clause
     * @param clause The vector of literals to add as a clause
     * 
     * @pre current decision level must be zero
     */
    void addExtDefClause(std::vector<CRef>& db, Lit ext_lit, const std::vector<Lit>& clause);

    /**
     * @brief Prioritize branching on a given set of variables
     * 
     * @param defs Extension variable definitions -- the extension variables will be prioritized
     */
    void prioritize(const std::vector<ExtDef>& defs);

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // EXTENSION VARIABLE SUBSTITUTION

    /**
     * @brief Check whether the given clause meets some condition and substitute extension
     * variables into a clause. Propagates substituted variables if they are unassigned.
     * 
     * @param clause The vector of literals in which to substitute
     * @param predicate The condition with which to check the clause
     */
  bool substitute(vec<Lit>& clause, SubstitutionPredicate& p);

    /**
     * @brief Check whether the given clause meets some condition and substitute extension
     * variables into a clause. Propagates substituted variables if they are unassigned.
     * 
     * @param clause The vector of literals in which to substitute
     */
    bool substitute(vec<Lit>& clause);

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // EXTENSION VARIABLE DELETION

    /**
     * @brief Get a set of extension variables to delete
     * 
     * @param varsToDelete The output set.
     * @param deletionPredicate A method for determining whether an extension variable should be removed.
     * 
     * @note Assumes that @code{varsToDelete} is initially empty
     */
    void getExtVarsToDelete(VarSet& varsToDelete, DeletionPredicate& deletionPredicate) const;

    /**
     * @brief Remove extension variables from the solver
     * 
     * @param setup a method to set up the data structures for @code{deletionPredicate}
     * @param deletionPredicate a method for determining whether an extension variable should be removed.
     */
    void deleteExtVars(DeletionPredicateSetup& setup, DeletionPredicate& deletionPredicate);

    /**
     * @brief Checks whether to delete variables and then calls @code{deleteExtVars}.
     * 
     * @param conflicts the current total number of conflicts seen by the solver
     */
    inline void checkDeleteExtVars(uint64_t conflicts);

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // DIP
  std::pair<bool,Lit> addDIPExtensionVariable(Lit a, Lit b);
  
  
public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // BUILT-IN HEURISTICS

    /**
     * @brief Add a variable according to the Local Extended Resolution heuristic
     * 
     */
    void generateLER(void);

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // STATISTICS

    /**
     * @brief Read the value of the extension timers
     * 
     * @param i 0: sel, 1: add, 2: delC, 3: delV, 4: sub, 5: stat
     * @return the value of the specified extension timer
     */
    double extTimerRead(unsigned int i) const;

    /**
     * @brief Start the overhead timer
     * 
     */
    void extTimerStart() const;

    /**
     * @brief Stop the overhead timer and add the time to the total overhead
     * 
     * @param ext_overhead The total overhead to add to
     */
    void extTimerStop(struct rusage& ext_overhead) const;

private:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // USER-DEFINED HEURISTICS

    // Decide whether to consider a clause based on its width
    bool user_extFilPredicate_width(CRef cr);
    // Decide whether to consider a clause based on its LBD
    bool user_extFilPredicate_lbd(CRef cr);
    // Decide whether to consider a clause based on the LER proof system (GlucosER)
    bool user_extFilPredicate_ler(CRef cr);
    
    // Select all clauses
    void user_extSelHeuristic_all     (std::vector<CRef>& output, const std::vector<CRef>& input, unsigned int numClauses);
    // Select clauses with highest activities
    void user_extSelHeuristic_activity(std::vector<CRef>& output, const std::vector<CRef>& input, unsigned int numClauses);

    // Define extension variables by selecting two literals at random
    void user_extDefHeuristic_random       (std::vector<ExtDef>& extVarDefBuffer, const std::vector<CRef>& selectedClauses, unsigned int maxNumNewVars);
    // Define extension variables by selecting the most frequently-appearing pairs of literals
    void user_extDefHeuristic_subexpression(std::vector<ExtDef>& extVarDefBuffer, const std::vector<CRef>& selectedClauses, unsigned int maxNumNewVars);
    // Define extension variables based on the LER proof system (GlucosER)
    void user_extDefHeuristic_ler          (std::vector<ExtDef>& extVarDefBuffer, const std::vector<CRef>& selectedClauses, unsigned int maxNumNewVars);
    
    // Decide whether to substitute into a clause based on its width and LBD
    bool user_extSubPredicate_size_lbd(vec<Lit>& clause);

    // Never delete extension variables
    void user_extDelPredicateSetup_none();
    bool user_extDelPredicate_none(Var x);
    // Always delete extension variables
    bool user_extDelPredicate_all(Var x);
    // Only delete extension variables with low activities
    void user_extDelPredicateSetup_activity();
    bool user_extDelPredicate_activity(Var x);

protected:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // HELPER FUNCTIONS FOR SUBSTITUTION

    /**
     * @brief Find the asserting definition clause for an extension variable
     * 
     * @param i_undef Return value - the index of the undefined literal in the clause
     * @param i_max Return value - the index of the literal in the clause with the highest decision level
     * @param x The unassigned extension literal
     * @param cs The list of potential asserting clauses - usually the list of definition clauses
     * @return The CRef for the asserting clause
     * 
     * @note pre-condition: the clause must be asserting at the current decision level
     */
    CRef findAssertingClause(int& i_undef, int& i_max, Lit x, std::vector<CRef>& cs);

protected:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // HELPER FUNCTIONS FOR CLAUSE DELETION

    /**
     * @brief Delete clauses from a learnt clause database if they contain variables that are being deleted
     * 
     * @param db the learnt clause database to delete from
     * @param varsToDelete the set of variables to delete (represented as a set of positive literals)
     */
    void deleteExtVarsFrom(vec<CRef>& db, VarSet& varsToDelete);

    /**
     * @brief Remove deleted clauses from CRef buffers. 
     * Used in tandem with @code{handleEventClauseDeleted} 
     * 
     * @note clears @code{m_deletedClauses}
     */
    inline void remove_flush();

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // EVENT HANDLERS

    /**
     * @brief Mark that a clause has been deleted -- used in tandem with @code{remove_flush}
     * 
     * @param cr the CRef of the deleted clause
     *
     * @note adds the CRef to @code{m_deletedClauses}
     */
    inline void handleEventClauseDeleted(CRef cr);

private:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // HELPER FUNCTIONS FOR USER-DEFINED HEURISTICS
    
    /**
     * @brief Find the literals that are different between two clauses
     * 
     * @param output return value -- a vector of literals that differ between the clauses
     * @param c1 the first clause
     * @param c2 the second clause
     * @return the number of literals that are different
     */
    int numDiffs(vec<Lit>& output, const Clause& c1, const Clause& c2);

    void quickselect_count(std::vector<LitPair>& db, std::tr1::unordered_map<LitPair, int>& subexpr_count, int l, int r, int k);
    std::tr1::unordered_map<LitPair, int> countSubexprs(std::vector<LitSet>& sets);
    std::vector<LitPair> getFreqSubexprs(std::tr1::unordered_map<LitPair, int>& subexpr_counts, unsigned int numSubexprs);
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION OF INLINE METHODS

///////////////////////
// STATE MODIFICATION

inline void ERManager::init(void) {
    originalNumVars = assignmentTrail.nVars();
}

inline void ERManager::newVar(Var v) {
    extensionLevel.push(0);
}

//////////////
// ACCESSORS

inline bool ERManager::isExtVar(Var x) const {
    return x >= originalNumVars;
}

inline bool ERManager::isCurrentExtVar(Var x) const {
    return x >= originalNumVars && xdm.containsExt(mkLit(x));
}

// FIXME: need to handle the case where the pair is currently queued for variable introduction in
// the extVarDefBuffer. We cannot just add directly to xdm since it might result in substitution
// before we introduce the new var.
inline bool ERManager::isValidDefPair(Lit a, Lit b, const std::tr1::unordered_set< std::pair<Lit, Lit> >& generatedPairs) const {
    // Ensure literal pair consists of different variables
    if (var(a) == var(b)) return false;

    // Ensure literals in pair are not set at level 0
    // FIXME: this breaks if the pair is over extension variables which have not been added yet
    if (assignmentTrail.value(a) != l_Undef && assignmentTrail.level(var(a)) == 0) return false;
    if (assignmentTrail.value(b) != l_Undef && assignmentTrail.level(var(b)) == 0) return false;
    
    // Ensure literal pair has not already been added
    if (xdm.containsPair(a, b)) return false;
    return generatedPairs.find(mkLitPair(a, b)) == generatedPairs.end();
}

///////////////
// STATISTICS

inline void ERManager::extTimerStart() const {
    getrusage(RUSAGE_SELF, &ext_timer_start);
}

inline void ERManager::extTimerStop(struct rusage& ext_overhead) const {
    getrusage(RUSAGE_SELF, &ext_timer_end);
    
    // Add to total overhead
    ext_overhead.ru_utime.tv_sec  += ext_timer_end.ru_utime.tv_sec - ext_timer_start.ru_utime.tv_sec;
    ext_overhead.ru_utime.tv_usec += ext_timer_end.ru_utime.tv_usec;

    // Check if subtracting the initial time would result in underflow
    if (ext_timer_start.ru_utime.tv_usec > ext_overhead.ru_utime.tv_usec) {
        ext_overhead.ru_utime.tv_usec += 1000000;
        ext_overhead.ru_utime.tv_sec  -= 1;
    }
    ext_overhead.ru_utime.tv_usec -= ext_timer_start.ru_utime.tv_usec;

    // Check if we carry over to the next second
    if (ext_overhead.ru_utime.tv_usec >= 1000000) {
        ext_overhead.ru_utime.tv_usec -= 1000000;
        ext_overhead.ru_utime.tv_sec  += 1;
    }
}

static inline double readTimer(const struct rusage& ext_overhead) {
    return (double)ext_overhead.ru_utime.tv_sec + (double)ext_overhead.ru_utime.tv_usec / 1000000;
}

inline double ERManager::extTimerRead(unsigned int i) const {
    switch(i) {
        case 0: return readTimer(ext_sel_overhead);
        case 1: return readTimer(ext_add_overhead);
        case 2: return readTimer(ext_delC_overhead);
        case 3: return readTimer(ext_delV_overhead);
        case 4: return readTimer(ext_sub_overhead);
        case 5: return readTimer(ext_stat_overhead);
        default: return -1.;
    }
}

/////////////////////
// CLAUSE SELECTION
inline void ERManager::filterIncremental(const CRef candidate, HeuristicType heuristicType) {
#if ER_ENABLE_GLUCOSER
    using namespace std::placeholders;
    static FilterPredicate fil_ler = std::bind(&ERManager::user_extFilPredicate_ler, this, _1);
    filterIncremental(candidate, fil_ler, HeuristicType::LER);
#endif
#if ER_USER_FILTER_HEURISTIC != ER_FILTER_HEURISTIC_NONE
    filterIncremental(candidate, user_extFilPredicate, heuristicType);
#endif
}

////////////////////////////////////
// EXTENSION VARIABLE INTRODUCTION

inline void ERManager::introduceExtVars(HeuristicType heuristicType) {
    introduceExtVars(extDefs, heuristicType);
}

inline void ERManager::addExtDefClause(std::vector<CRef>& db, Lit ext_lit, const std::initializer_list<Lit>& clause) {
    tmp_vec.clear(); for (const auto l : clause) tmp_vec.push(l);
    addExtDefClause(db, ext_lit, tmp_vec);
}

inline void ERManager::addExtDefClause(std::vector<CRef>& db, Lit ext_lit, const std::vector<Lit>& clause) {
    tmp_vec.clear(); for (const auto l : clause) tmp_vec.push(l);
    addExtDefClause(db, ext_lit, tmp_vec);
}

inline void ERManager::checkGenerateDefinitions(uint64_t conflicts) {
    if (conflicts - prevExtensionConflict >= static_cast<unsigned int>(ext_freq)) {
        prevExtensionConflict = conflicts;
        selectClauses(user_extSelHeuristic, ERManager::HeuristicType::DEFAULT);
        defineExtVars(user_extDefHeuristic);
    }
}

////////////////////////////////////
// EXTENSION VARIABLE SUBSTITUTION

/**
 * @brief Check whether the given clause meets some condition and substitute extension
 * variables into a clause. Propagates substituted variables if they are unassigned.
 * 
 * @param clause The vector of literals in which to substitute
 */
  inline bool ERManager::substitute(vec<Lit>& clause) {
#if ER_USER_SUBSTITUTE_HEURISTIC != ER_SUBSTITUTE_HEURISTIC_NONE
    return substitute(clause, user_extSubPredicate);
#endif
    return false;
}

////////////////////////////////
// EXTENSION VARIABLE DELETION

inline void ERManager::checkDeleteExtVars(uint64_t conflicts) {
    if (conflicts - prevDelExtVarConflict >= static_cast<unsigned int>(ext_del_freq)){
        prevDelExtVarConflict = conflicts;
        deleteExtVars(user_extDelPredicateSetup, user_extDelPredicate);
    }
}

inline static void removeSet(std::vector<CRef>& vec, const std::tr1::unordered_set<CRef>& toDelete) {
    unsigned int i, j;
    for (i = j = 0; i < vec.size(); i++)
        if (toDelete.find(vec[i]) == toDelete.end())
            vec[j++] = vec[i];
    vec.erase(vec.begin() + j, vec.end());
}

inline void ERManager::remove_flush() {
    // Remove CRefs from buffers of filtered clauses
    removeSet(m_filteredClauses, m_deletedClauses);
    removeSet(m_filteredClauses_ler, m_deletedClauses);

    // Remove CRefs from buffer of selected clauses
    removeSet(m_selectedClauses, m_deletedClauses);
    removeSet(m_selectedClauses_ler, m_deletedClauses);

    // Clear deletion buffer
    m_deletedClauses.clear();
}

///////////////////
// EVENT HANDLERS

inline void ERManager::handleEventClauseDeleted(CRef cr) {
    m_deletedClauses.insert(cr);
}

////////////////////////
// BUILT-IN HEURISTICS

inline void ERManager::generateLER(void) {
    using namespace std::placeholders;
    static SelectionHeuristic sel_ler = std::bind(&ERManager::user_extSelHeuristic_all, this, _1, _2, _3);
    static ExtDefHeuristic def_ler = std::bind(&ERManager::user_extDefHeuristic_ler, this, _1, _2, _3);

    // Generate extension variable definitions
    selectClauses(sel_ler, HeuristicType::LER, 1);
    defineExtVars(def_ler, HeuristicType::LER);
}

}

#endif
