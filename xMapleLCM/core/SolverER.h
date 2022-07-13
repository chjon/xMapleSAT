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

    inline bool isExtVar(Var x) const;

#ifdef TESTING
    inline void set_value(Var x, lbool v, int l);
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

    // Extension Variable Introduction
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
    void enforceWatcherInvariant(vec<Lit>& clause);

    // Extension Variable Substitution

    /**
     * @brief Check whether the given clause meets some condition and substitute extension variables into a clause
     * 
     * @param clause The vector of literals in which to substitute
     * @param predicate The condition with which to check the clause
     */
    inline void substitute(vec<Lit>& clause, SubstitutionPredicate& p) const;

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

    FilterPredicate       user_extFilPredicate;
    SelectionHeuristic    user_extSelHeuristic;
    ExtDefHeuristic       user_extDefHeuristic;
    SubstitutionPredicate user_extSubPredicate;

    ////////////////
    // Statistics //
    ////////////////

    uint64_t total_ext_vars, deleted_ext_vars, max_ext_vars;
    uint64_t conflict_extclauses, learnt_extclauses, lbd_total, branchOnExt;

    double extTimerRead(unsigned int i); // 0: sel, 1: add, 2: delC, 3: delV, 4: sub, 5: stat

protected:

    // // Update stats
    // void updateExtFracStat(vec<Lit>& clause) {
    //     int numExtVarsInClause = getNumExtVars(clause);
    //     double extFrac = numExtVarsInClause / (double) clause.size();
    //     extfrac_total += extFrac;
    // }

    // Measuring extended resolution overhead
    struct rusage ext_timer_start, ext_timer_end;
    struct rusage ext_sel_overhead; // Overhead for selecting clauses for adding extension variables
    struct rusage ext_add_overhead; // Overhead for adding extension variables
    struct rusage ext_delC_overhead; // Overhead for deleting clauses containing extension variables
    struct rusage ext_delV_overhead; // Overhead for deleting extension variables
    struct rusage ext_sub_overhead; // Overhead for substituting disjunctions containing extension variables
    struct rusage ext_stat_overhead; // Overhead for measuring statistics

    void   extTimerStart();
    void   extTimerStop(struct rusage& ext_overhead);

    Solver* solver;
    ExtDefMap<Lit> xdm;

    std::vector<CRef> m_filteredClauses;
    std::vector<CRef> m_selectedClauses;
    std::vector<ExtDef> m_extVarDefBuffer;
    vec<Lit> tmp;

#ifdef TESTING
    std::map< Var, std::pair<lbool, int> > test_value;
#endif

};

#ifdef TESTING
void  SolverER::set_value(Var x, lbool v, int l) { test_value.insert(std::make_pair(x, std::make_pair(v, l))); }
int   SolverER::level(Var x) const { auto it = test_value.find(x); return (it == test_value.end()) ? (-1     ) : (it->second.second); }
lbool SolverER::value(Var x) const { auto it = test_value.find(x); return (it == test_value.end()) ? (l_Undef) : (it->second.first ); }
lbool SolverER::value(Lit p) const { return value(var(p)) ^ sign(p); }
#else
int   SolverER::level(Var x) const { return solver->level(x); }
lbool SolverER::value(Var x) const { return solver->value(x); }
lbool SolverER::value(Lit p) const { return solver->value(p); }
#endif

void SolverER::substitute(vec<Lit>& clause, SubstitutionPredicate& p) const {
    if (p(clause)) xdm.substitute(clause);
}

// EXTENDED RESOLUTION - statistics
inline void SolverER::extTimerStart() {
    getrusage(RUSAGE_SELF, &ext_timer_start);
}

inline void SolverER::extTimerStop(struct rusage& ext_overhead) {
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
    return x > originalNumVars;
}

}

#endif