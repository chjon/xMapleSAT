/*******************************************************************************************[SolverER.h]
Copyright (c) 2022, Jonathan Chung

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

#ifndef Minisat_SolverER_h
#define Minisat_SolverER_h

#include <sys/resource.h>
#include <functional>
#include <initializer_list>
#include <map>
#include <tuple>
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#include <utility>
#include <vector>
#include <mtl/Vec.h>
#include <mtl/ExtDefMap.h>
#include <core/Solver.h>
#include <core/SolverTypes.h>
#include <core/SolverERTypes.h>

namespace Minisat {

class SolverER {
public:
    SolverER(Solver* s);
    ~SolverER();

    inline int   level(Var x) const;
    inline lbool value(Var x) const;
    inline lbool value(Lit p) const;

    int originalNumVars; // The number of variables in the original formula
                         // This value is used to quickly check whether a variable is an extension variable
    std::tr1::unordered_map<Var, std::vector<CRef> > extDefs; // List of extension definition clauses.

    // Determine whether a variable is an extension variable
    inline bool isExtVar(Var x) const;
    
    // Determine whether a variable is an extension variable in the data structure for variable substitution
    inline bool isCurrentExtVar(Var x) const;

    // Determine whether a pair of literals can be used as the basis literals for a new extension variable
    inline bool isValidDefPair(Lit a, Lit b, const std::tr1::unordered_set<LitPair>& generatedPairs) const;

#ifdef TESTING
    inline void set_value(Var x, lbool v, int l);
    inline void addToExtDefBuffer(ExtDef extDef);
    inline unsigned int extDefBufferSize();
#endif

    // Clause Selection

    /**
     * @brief Filter clauses for clause selection
     * 
     * @param candidates a list of clauses to check
     * @param filterPredicate a method for selecting clauses for defining extension variables
     * 
     * @note Saves outputs in @code{m_filteredClauses}
     */
    void filterBatch(const vec<CRef>& candidates, FilterPredicate& filterPredicate);

    /**
     * @brief Filter clauses for clause selection
     * 
     * @param candidate the clause to check
     * @param filterPredicate a method for selecting clauses for defining extension variables
     * 
     * @note Saves outputs in @code{m_filteredClauses}
     */
    void filterIncremental(const CRef candidate, FilterPredicate& filterPredicate);

    /**
     * @brief Select clauses for defining extension variables
     * 
     * @param selectionHeuristic a method for selecting clauses for defining extension variables
     * 
     * @note Takes input from @code{m_filteredClauses}
     * @note Saves outputs in @code{m_selectedClauses}
     */
    void selectClauses(SelectionHeuristic& selectionHeuristic);

    // Extension Variable Definition

    /**
     * @brief Generate extension variable definitions
     * 
     * @param extDefHeuristic a method for generating extension variable definitions given a list of selected clauses
     */
    void defineExtVars(ExtDefHeuristic& extDefHeuristic);
    
    /**
     * @brief Checks whether to generate definitions and then calls @code{selectClauses} and @code{defineExtVars}.
     */
    inline void generateDefinitions();

    // Extension Variable Introduction

    /**
     * @brief Introduce extension variables into the solver
     * 
     * @param ext_def_db The map of extension variables to a list of their corresponding extension definition clauses.
     * In most cases this should be equal to @code{extDefs}.
     */
    void introduceExtVars(std::tr1::unordered_map<Var, std::vector<CRef> >& ext_def_db);

    /**
     * @brief Adds an extension definition clause to the appropriate clause database
     * 
     * @param db The clause database to which to add the clause
     * @param ext_lit The extension variable corresponding to the given clause
     * @param clause The vector of literals to add as a clause
     * 
     * @note Condition: current decision level must be zero
     */
    void addExtDefClause(std::vector<CRef>& db, Lit ext_lit, vec<Lit>& clause);
    void addExtDefClause(std::vector<CRef>& db, Lit ext_lit, const std::initializer_list<Lit>& clause);
    void addExtDefClause(std::vector<CRef>& db, Lit ext_lit, const std::vector<Lit>& clause);

    /**
     * @brief Prioritize branching on a given set of variables
     * 
     * @param defs Extension variable definitions -- the extension variables will be prioritized
     */
    void prioritize(const std::vector<ExtDef>& defs);

    /**
     * @brief Ensures the first two literals are in the right order for the watchers
     * 
     * @param clause a vector of multiple literals
     * 
     * @note @code{clause} must have length > 1. Unary clauses should be propagated directly and not learnt
     */
    void enforceWatcherInvariant(vec<Lit>& clause) const;

    // Extension Variable Substitution

    /**
     * @brief Check whether the given clause meets some condition and substitute extension variables into a clause
     * 
     * @param clause The vector of literals in which to substitute
     * @param predicate The condition with which to check the clause
     * @return true if a variable was substituted into the clase
     * @return false otherwise
     */
    inline bool substitute(vec<Lit>& clause, SubstitutionPredicate& predicate) const;

    /**
     * @brief Enforce the invariant that learnt clauses are asserting after backtracking. Rarely, after variable
     * substitution and backtracking, some extension variables are undefined even though their definitions are
     * falsified.
     * 
     * @param clause The learnt clause to check
     */
    void enforceLearntClauseInvariant(const vec<Lit>& clause);

    // Extension Variable Deletion

    /**
     * @brief Get a set of extension variables to delete
     * 
     * @param varsToDelete The output set.
     * @param deletionPredicate A method for determining whether an extension variable should be removed.
     * 
     * @note Assumes that @code{varsToDelete} is initially empty
     */
    void getExtVarsToDelete(LitSet& varsToDelete, DeletionPredicate& deletionPredicate) const;

    /**
     * @brief Remove extension variables from the solver
     * 
     * @param deletionPredicate a method for determining whether an extension variable should be removed.
     */
    void deleteExtVars(DeletionPredicate& deletionPredicate);
    
    ///////////////////////
    // Solver.h helpers //
    ///////////////////////
    // Since SolverER has its own data structures, we implement these methods here for readability in the main solver code

    /**
     * @brief Relocate CRefs to new ClauseAllocator
     * 
     * @param to The ClauseAllocator into which to reloc 
     */
    void relocAll(ClauseAllocator& to);

    // TODO: verify safety
    void removeSatisfied();

    /////////////////////////////
    // USER-DEFINED HEURISTICS //
    /////////////////////////////
    bool user_extFilPredicate_width(CRef cr);
    bool user_extFilPredicate_lbd(CRef cr);
    
    // Select all clauses
    void user_extSelHeuristic_all     (std::vector<CRef>& output, const std::vector<CRef>& input, unsigned int numClauses);
    // Select clauses with highest activities
    void user_extSelHeuristic_activity(std::vector<CRef>& output, const std::vector<CRef>& input, unsigned int numClauses);
    
    void user_extDefHeuristic_random       (std::vector<ExtDef>& extVarDefBuffer, const std::vector<CRef>& selectedClauses, unsigned int maxNumNewVars);
    void user_extDefHeuristic_subexpression(std::vector<ExtDef>& extVarDefBuffer, const std::vector<CRef>& selectedClauses, unsigned int maxNumNewVars);
    
    bool user_extSubPredicate_size_lbd(vec<Lit>& clause);

    // Never delete extension variables
    bool user_extDelPredicate_none(Var x);
    // Always delete extension variables
    bool user_extDelPredicate_all(Var x);
    // Only delete extension variables with low activities
    bool user_extDelPredicate_activity(Var x);

    FilterPredicate       user_extFilPredicate;
    SelectionHeuristic    user_extSelHeuristic;
    ExtDefHeuristic       user_extDefHeuristic;
    SubstitutionPredicate user_extSubPredicate;
    DeletionPredicate     user_extDelPredicate;

    ////////////////
    // Statistics //
    ////////////////

    uint64_t total_ext_vars, deleted_ext_vars, max_ext_vars;
    uint64_t conflict_extclauses, learnt_extclauses, lbd_total, branchOnExt;

    double extTimerRead(unsigned int i); // 0: sel, 1: add, 2: delC, 3: delV, 4: sub, 5: stat

protected:
    long unsigned int prevExtensionConflict; // Stores the last time extension variables were added

    /////////////////////////////
    // Command-line parameters //
    /////////////////////////////
    int       ext_freq;           // Number of conflicts to wait before trying to introduce an extension variable              (default 2000)
    int       ext_window;         // Number of clauses to consider when introducing extension variables.                       (default 100)
    int       ext_max_intro;      // Maximum number of extension variables to introduce at once.                               (default 1)
    double    ext_prio_act;       // The fraction of maximum activity that should be given to new variables                    (default 0.5)
    bool      ext_pref_sign;      // Preferred sign for new variables                                                          (default true (negated))
    int       ext_min_width;      // Minimum clause width to consider when selecting clauses
    int       ext_max_width;      // Maximum clause width to consider when selecting clauses
    int       ext_filter_num;     // Maximum number of clauses after the filter step
    int       ext_sub_min_width;  // Minimum width of clauses to substitute into
    int       ext_sub_max_width;  // Maximum width of clauses to substitute into
    int       ext_min_lbd;        // Minimum LBD of clauses to substitute into
    int       ext_max_lbd;        // Maximum LBD of clauses to substitute into
    double    ext_act_threshold;  // Activity threshold for deleting clauses

    // // Update stats
    // void updateExtFracStat(vec<Lit>& clause) {
    //     int numExtVarsInClause = getNumExtVars(clause);
    //     double extFrac = numExtVarsInClause / (double) clause.size();
    //     extfrac_total += extFrac;
    // }

    // Measuring extended resolution overhead
    mutable struct rusage ext_timer_start, ext_timer_end;
    mutable struct rusage ext_sel_overhead; // Overhead for selecting clauses for adding extension variables
    mutable struct rusage ext_add_overhead; // Overhead for adding extension variables
    mutable struct rusage ext_delC_overhead; // Overhead for deleting clauses containing extension variables
    mutable struct rusage ext_delV_overhead; // Overhead for deleting extension variables
    mutable struct rusage ext_sub_overhead; // Overhead for substituting disjunctions containing extension variables
    mutable struct rusage ext_stat_overhead; // Overhead for measuring statistics

    void   extTimerStart() const;
    void   extTimerStop(struct rusage& ext_overhead) const;

    Solver* solver;
    ExtDefMap<Lit> xdm;

    std::vector<CRef> m_filteredClauses;
    std::vector<CRef> m_selectedClauses;
    std::vector<ExtDef> m_extVarDefBuffer;

    // Keep track of clauses which have been removed so that they can be removed from the buffers above
    std::tr1::unordered_set<CRef> m_deletedClauses;
    vec<Lit> tmp;

    /////////////////////////////////
    // HELPERS FOR CLAUSE DELETION //
    /////////////////////////////////
    inline void remove_incremental(CRef cr);
    inline void remove_flush();

#ifdef TESTING
    std::map< Var, std::pair<lbool, int> > test_value;
#endif

private:
    /////////////////////////////////////////
    // HELPERS FOR USER-DEFINED HEURISTICS //
    /////////////////////////////////////////
    void quickselect_count(std::vector<LitPair>& db, std::tr1::unordered_map<LitPair, int>& subexpr_count, int l, int r, int k);
    std::tr1::unordered_map<LitPair, int> countSubexprs(std::vector<LitSet>& sets);
    std::vector<LitPair> getFreqSubexprs(std::tr1::unordered_map<LitPair, int>& subexpr_counts, unsigned int numSubexprs);
};

#ifdef TESTING
void  SolverER::set_value(Var x, lbool v, int l) {
    auto it = test_value.find(x);
    if (it == test_value.end()) test_value.insert(std::make_pair(x, std::make_pair(v, l)));
    else it->second = std::make_pair(v, l);
}
void  SolverER::addToExtDefBuffer(ExtDef extDef) { m_extVarDefBuffer.push_back(extDef); }
unsigned int SolverER::extDefBufferSize() { return m_extVarDefBuffer.size(); }

int   SolverER::level(Var x) const { auto it = test_value.find(x); return (it == test_value.end()) ? (-1     ) : (it->second.second); }
lbool SolverER::value(Var x) const { auto it = test_value.find(x); return (it == test_value.end()) ? (l_Undef) : (it->second.first ); }
lbool SolverER::value(Lit p) const { return value(var(p)) ^ sign(p); }
#else
int   SolverER::level(Var x) const { return solver->level(x); }
lbool SolverER::value(Var x) const { return solver->value(x); }
lbool SolverER::value(Lit p) const { return solver->value(p); }
#endif

bool SolverER::substitute(vec<Lit>& clause, SubstitutionPredicate& p) const {
    bool subbed = false;
    extTimerStart();
    // xdm.absorb(clause);
    if (p(clause)) {
        vec<Lit> extLits;
        xdm.substitute(clause, extLits);
        subbed = extLits.size() > 0;
    }
    extTimerStop(ext_sub_overhead);
    return subbed;
}

// EXTENDED RESOLUTION - statistics
inline void SolverER::extTimerStart() const {
    getrusage(RUSAGE_SELF, &ext_timer_start);
}

inline void SolverER::extTimerStop(struct rusage& ext_overhead) const {
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

static inline double readTimer(struct rusage& ext_overhead) {
    return (double)ext_overhead.ru_utime.tv_sec + (double)ext_overhead.ru_utime.tv_usec / 1000000;
}
inline double SolverER::extTimerRead(unsigned int i) {
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

inline void SolverER::addExtDefClause(std::vector<CRef>& db, Lit ext_lit, const std::initializer_list<Lit>& clause) {
    tmp.clear(); for (const auto l : clause) tmp.push(l);
    addExtDefClause(db, ext_lit, tmp);
}
inline void SolverER::addExtDefClause(std::vector<CRef>& db, Lit ext_lit, const std::vector<Lit>& clause) {
    tmp.clear(); for (const auto l : clause) tmp.push(l);
    addExtDefClause(db, ext_lit, tmp);
}

inline bool SolverER::isExtVar(Var x) const {
    return x >= originalNumVars;
}

inline bool SolverER::isCurrentExtVar(Var x) const {
    return x >= originalNumVars && xdm.containsExt(mkLit(x));
}

// FIXME: need to handle the case where the pair is currently queued for variable introduction in the extVarDefBuffer
// We cannot just add directly to xdm since it might result in substitution before we introduce the new var
inline bool SolverER::isValidDefPair(Lit a, Lit b, const std::tr1::unordered_set< std::pair<Lit, Lit> >& generatedPairs) const {
    // Ensure literal pair consists of different variables
    if (var(a) == var(b)) return false;

    // Ensure literals in pair are not set at level 0
    if (value(a) != l_Undef && level(var(a)) == 0) return false;
    if (value(b) != l_Undef && level(var(b)) == 0) return false;
    
    // Ensure literal pair has not already been added
    if (xdm.containsPair(a, b)) return false;
    return generatedPairs.find(mkLitPair(a, b)) == generatedPairs.end();
}

inline void SolverER::generateDefinitions() {
    if (solver->conflicts - prevExtensionConflict >= static_cast<unsigned int>(ext_freq)) {
        prevExtensionConflict = solver->conflicts;
        selectClauses(user_extSelHeuristic);
        defineExtVars(user_extDefHeuristic);
    }
}

inline void SolverER::remove_incremental(CRef cr) {
    m_deletedClauses.insert(cr);
}

inline void SolverER::remove_flush() {
    unsigned int i, j;

    // Remove CRefs from buffer of filtered clauses
    for (i = j = 0; i < m_filteredClauses.size(); i++) {
        if (m_deletedClauses.find(m_filteredClauses[i]) == m_deletedClauses.end()) {
            m_filteredClauses[j++] = m_filteredClauses[i];
        }
    }
    m_filteredClauses.erase(m_filteredClauses.begin() + j, m_filteredClauses.end());

    // Remove CRefs from buffer of selected clauses
    for (i = j = 0; i < m_selectedClauses.size(); i++) {
        if (m_deletedClauses.find(m_selectedClauses[i]) == m_deletedClauses.end()) {
            m_selectedClauses[j++] = m_selectedClauses[i];
        }
    }
    m_selectedClauses.erase(m_selectedClauses.begin() + j, m_selectedClauses.end());

    // Clear deletion buffer
    m_deletedClauses.clear();
}

}

#endif