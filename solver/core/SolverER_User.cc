#include <algorithm>
#include "core/Solver.h"

namespace Minisat {

static inline std::pair<Lit, Lit> mkLitPair(Lit a, Lit b) {
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

struct clause_width_lt { 
    ClauseAllocator& ca;
    clause_width_lt(ClauseAllocator& ca_) : ca(ca_) {}
    bool operator () (CRef x, CRef y) { 
        return ca[x].size() < ca[y].size();
    }
};

void Solver::user_er_filter_incremental(const CRef candidate) {
    // Filter clauses based on their sizes
#if ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_RANGE
    const int k = ca[candidate].size();
    if (k >= ext_min_width && k <= ext_max_width) er_filteredClauses.push_back(candidate);
#elif ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LONGEST
    struct clause_width_lt lt = clause_width_lt(ca);
    if (er_filteredClauses.size() < static_cast<unsigned int>(ext_filter_num)) {
        er_filteredClauses.push_back(candidate);
        std::push_heap(er_filteredClauses.begin(), er_filteredClauses.end(), lt);
    } else {
        const int k = ca[candidate].size();
        int old_min = ca[er_filteredClauses.back()].size();
        if (k > old_min) {
            er_filteredClauses.push_back(candidate);
            std::push_heap(er_filteredClauses.begin(), er_filteredClauses.end(), lt);
            std::pop_heap (er_filteredClauses.begin(), er_filteredClauses.end(), lt);
        }
    }
#elif ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LBD
    // Filter clauses based on their creation LBD
    if (ca[candidate].good_lbd()) er_filteredClauses.push_back(candidate);
#elif ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_GLUCOSER
    if (ca[candidate].learnt()) {
        if (er_filteredClauses.size() > 1) {
            er_filteredClauses[0] = er_filteredClauses[1];
            er_filteredClauses[1] = candidate;
        } else {
            er_filteredClauses.push_back(candidate);
        }
    }
#endif
}

void Solver::user_er_filter_delete_incremental(CRef cr) {
    extTimerStart();
#if ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_RANGE   || \
    ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LBD     || \
    ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LONGEST || \
    ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_GLUCOSER
    er_deletedClauses.insert(cr);
#endif
    extTimerStop(ext_delC_overhead);
}

void Solver::user_er_filter_delete_flush(void) {
#if ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_RANGE || \
    ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LBD     || \
    ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LONGEST || \
    ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_GLUCOSER
    extTimerStart();
    std::vector<CRef> tmp;
    for (std::vector<CRef>::iterator it = er_filteredClauses.begin(); it != er_filteredClauses.end(); it++)
        if (er_deletedClauses.find(*it) == er_deletedClauses.end())
            tmp.push_back(*it);
    er_filteredClauses = tmp;
    er_deletedClauses.clear();
#if ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LONGEST
    std::make_heap(er_filteredClauses);
#endif
    extTimerStop(ext_delC_overhead);
#endif
}

void Solver::user_er_reloc(ClauseAllocator& to) {
#if ER_USER_FILTER_HEURISTIC != ER_FILTER_HEURISTIC_NONE
    extTimerStart();
    for (unsigned int i = 0; i < er_filteredClauses.size(); i++) {
        ca.reloc(er_filteredClauses[i], to);
    }
    extTimerStop(ext_sel_overhead);
#endif
}

std::vector<CRef> Solver::user_er_select(Solver& solver, unsigned int numClauses) {
#if ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_NONE
    return user_er_select_naive(solver, numClauses);
#elif ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_ACTIVITY
    return user_er_select_activity(solver, numClauses);
#elif ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_ACTIVITY2
    return user_er_select_activity2(solver, numClauses);
#elif ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_GLUCOSER
    return user_er_select_glucosER(solver, numClauses);
#endif
}

#if ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_NONE
std::vector<CRef> Solver::user_er_select_naive(Solver& s, unsigned int numClauses) {
    std::vector<CRef> clauseWindow;
#if ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_RANGE || ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LBD || ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LONGEST
    for (std::vector<CRef>::iterator it = s.er_filteredClauses.begin(); it != s.er_filteredClauses.end(); it++) {
        if (numClauses == 0) break;
        clauseWindow.push_back(*it);
        numClauses--;
    }
#endif
    return clauseWindow;
}

#elif ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_ACTIVITY
// Build a window of the clauses with the top k highest activities
static void addClauseToWindow(ClauseAllocator& ca, std::vector<CRef>& window, CRef clauseIndex, unsigned int maxWindowSize) {
    CRef tmp = 0;
    double clauseActivity = ca[clauseIndex].activity();
    for (unsigned int i = 0; i < window.size() && i < maxWindowSize; i++) {
        if (clauseIndex == window[i]) return;
        if (clauseActivity > ca[window[i]].activity()) {
            tmp = window[i];
            window[i] = clauseIndex;
            clauseActivity = ca[tmp].activity();
            clauseIndex = tmp;
        }
    }
    if (window.size() != maxWindowSize) {
        window.push_back(clauseIndex);
    }
}

// EXTENDED RESOLUTION - clause selection heuristic
std::vector<CRef> Solver::user_er_select_activity(Solver& s, unsigned int numClauses) {
    // Find the variables in the clauses with the top k highest activities
    std::vector<CRef> clauseWindow;

    // Use incremental filtered clause list
#if ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_RANGE || ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LBD
    for (std::vector<CRef>::iterator it = s.er_filteredClauses.begin(); it != s.er_filteredClauses.end(); it++)
        addClauseToWindow(s.ca, clauseWindow, *it, numClauses);
#endif

    return clauseWindow;
}

#elif ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_ACTIVITY2
// Partition elements such that clauses with larger activities are left of the pivot
// This is a helper function for quickselect
inline int Solver::partition_activity(vec<CRef>& db, ClauseAllocator& ca, int l, int r, int pivot) {
    const double pivotActivity = ca[db[pivot]].activity();
    CRef tmp = 0;
    tmp = db[pivot]; db[pivot] = db[r]; db[r] = tmp;
    int i = l;

    for (int j = l; j <= r - 1; j++) {
        if (ca[db[j]].activity() > pivotActivity) {
            // Swap db[i] and db[j]
            tmp = db[i]; db[i] = db[j]; db[j] = tmp;
            i++;
        }
    }
    // Swap db[i] and db[r]
    tmp = db[i]; db[i] = db[r]; db[r] = tmp;
    return i;
}

// Move elements with the k largest activities to the front of the list
void Solver::quickselect_activity(vec<CRef>& db, Solver& solver, int l, int r, int k) {
    // Ensure we have a valid value of k
    k--;
    if (k <= 0 || k >= r - l + 1) return;
    while (!solver.interrupted()) {
        // Partition the array around last element and get position of pivot element in sorted array
        int pivot = l + solver.irand(solver.random_seed, r - l + 1);
        pivot = partition_activity(db, solver.ca, l, r, pivot);

        // Update selection bounds
        if      (k < pivot) r = pivot - 1;
        else if (k > pivot) l = pivot + 1;
        else break;
    }
}

// Get the elements with the k largest activities and store them in the target vector
void Solver::copy_k_largest_activity(vec<CRef>& target, vec<CRef>& source, Solver& solver, unsigned int k) {
#define ER_ALLOW_DB_MODIFICATION false
#if ER_ALLOW_DB_MODIFICATION
    quickselect_activity(source, solver, 0, source.size() - 1, k);
    for (unsigned int i = 0; i < std::min(static_cast<unsigned int>(source.size()), k); i++) target.push(source[i]);
#else
    // TODO: is making this copy necessary?
    vec<CRef> copy;
    source.copyTo(copy);
    quickselect_activity(copy, solver, 0, copy.size() - 1, k);
    for (unsigned int i = 0; i < std::min(static_cast<unsigned int>(copy.size()), k); i++) target.push(copy[i]);
#endif
}

std::vector<CRef> Solver::user_er_select_activity2(Solver& s, unsigned int numClauses) {
    // Find the variables in the clauses with the top k highest activities
    // Uses a quicksort-type of algo for expected-linear time
    // Adapted from: https://en.wikipedia.org/wiki/Quickselect
    // If this still takes too long, consider: https://en.wikipedia.org/wiki/Floyd-Rivest_algorithm

    vec<CRef> clauses;
    vec<CRef> filteredClauses;

    // Use incremental filtered clause list
#if ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_RANGE || ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LBD || ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LONGEST
    for (std::vector<CRef>::iterator it = s.er_filteredClauses.begin(); it != s.er_filteredClauses.end(); it++)
        filteredClauses.push(*it);
#endif
    copy_k_largest_activity(clauses, filteredClauses, s, numClauses);

    std::vector<CRef> clauseWindow;
    quickselect_activity(clauses, s, 0, clauses.size() - 1, numClauses);
    for (unsigned int i = 0; i < std::min(static_cast<unsigned int>(clauses.size()), numClauses); i++) clauseWindow.push_back(clauses[i]);
    return clauseWindow;
}

#elif ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_GLUCOSER
std::vector<CRef> Solver::user_er_select_glucosER(Solver& s, unsigned int numClauses) {
    return s.er_filteredClauses;
}
#endif

// EXTENDED RESOLUTION - variable definition heuristic
#if ER_USER_ADD_HEURISTIC != ER_ADD_HEURISTIC_NONE
std::vector< std::pair< Var, std::pair<Lit, Lit> > > Solver::user_er_add(Solver& s, std::vector<CRef>& candidateClauses, unsigned int maxNumNewVars) {
#if ER_USER_ADD_HEURISTIC == ER_ADD_HEURISTIC_SUBEXPR
    return user_er_add_subexpr(s, candidateClauses, maxNumNewVars);
#elif ER_USER_ADD_HEURISTIC == ER_ADD_HEURISTIC_RANDOM
    return user_er_add_random(s, candidateClauses, maxNumNewVars);
#elif ER_USER_ADD_HEURISTIC == ER_ADD_HEURISTIC_GLUCOSER
    return user_er_add_glucosER(s, candidateClauses, maxNumNewVars);
#endif
}
#endif

#if ER_USER_ADD_HEURISTIC == ER_ADD_HEURISTIC_SUBEXPR
// Partition elements such that clauses with larger activities are left of the pivot
// This is a helper function for quickselect
inline int Solver::partition_count(std::vector< std::pair<Lit, Lit> >& db, std::tr1::unordered_map<std::pair<Lit, Lit>, int>& subexpr_count, int l, int r, int pivot) {
    const int pivot_count = subexpr_count.find(db[pivot])->second;
    std::pair<Lit, Lit> tmp;
    tmp = db[pivot]; db[pivot] = db[r]; db[r] = tmp;
    int i = l;

    for (int j = l; j <= r - 1; j++) {
        if (subexpr_count.find(db[j])->second > pivot_count) {
            // Swap db[i] and db[j]
            tmp = db[i]; db[i] = db[j]; db[j] = tmp;
            i++;
        }
    }
    // Swap db[i] and db[r]
    tmp = db[i]; db[i] = db[r]; db[r] = tmp;
    return i;
}

void Solver::quickselect_count(std::vector< std::pair<Lit, Lit> >& db, std::tr1::unordered_map<std::pair<Lit, Lit>, int>& subexpr_count, Solver& solver, int l, int r, int k) {
    // Ensure we have a valid value of k
    k--;
    if (k <= 0 || k >= r - l + 1) return;
    while (!solver.interrupted()) {
        // Partition the array around last element and get position of pivot element in sorted array
        int pivot = l + solver.irand(solver.random_seed, r - l + 1);
        pivot = partition_count(db, subexpr_count, l, r, pivot);

        // Update selection bounds
        if      (k < pivot) r = pivot - 1;
        else if (k > pivot) l = pivot + 1;
        else break;
    }
}

static inline void addIntersectionToSubexprs(std::tr1::unordered_map<std::pair<Lit, Lit>, int>& subexprs, const std::vector<Lit>& intersection) {
    // Time complexity: O(k^2)
    for (unsigned int i = 0; i < intersection.size(); i++) {
        for (unsigned int j = i + 1; j < intersection.size(); j++) {
            // Count subexpressions of length 2
            std::pair<Lit, Lit> key = mkLitPair(intersection[i], intersection[j]);

            // Add to the counter for this literal pair
            std::tr1::unordered_map<std::pair<Lit, Lit>, int>::iterator it = subexprs.find(key);
            if (it == subexprs.end()) subexprs.insert(std::make_pair(key, 1));
            else it->second++;
        }
    }
}

static inline std::tr1::unordered_map<std::pair<Lit, Lit>, int> countSubexprs(const Solver& s, std::vector< std::tr1::unordered_set<Lit> >& sets) {
    // Count subexpressions by looking at intersections
    // Time complexity: O(w k^2)
    std::tr1::unordered_map<std::pair<Lit, Lit>, int> subexprs;

#if ER_USER_ADD_SUBEXPR_SET_INTERSECTION
    for (unsigned int i = 0; i < sets.size(); i++) {
        for (unsigned int j = i + 1; j < sets.size(); j++) {
            // TODO: Check if we've already processed a pair of clauses (cache) and add their counts

            // We might spend a lot of time here - exit if interrupted
            // FIXME: ideally, we wouldn't have to check for this at all if the sets of literals were sufficiently small
            if (s.interrupted()) goto SUBEXPR_DOUBLE_BREAK;
            std::vector<Lit> intersection(sets[i].size() + sets[j].size());
            std::vector<Lit>::iterator it = std::set_intersection(sets[i].begin(), sets[i].end(), sets[j].begin(), sets[j].end(), intersection.begin());
            intersection.resize(it - intersection.begin());
            addIntersectionToSubexprs(subexprs, intersection);
        }
    }
#else
    for (unsigned int i = 0; i < sets.size(); i++) {
        std::tr1::unordered_set<Lit>& clause = sets[i];
        for (std::tr1::unordered_set<Lit>::iterator j = clause.begin(); j != clause.end(); j++) {
            std::tr1::unordered_set<Lit>::iterator k = j; k++;
            // We might spend a lot of time here - exit if interrupted
            // FIXME: ideally, we wouldn't have to check for this at all if the sets of literals were sufficiently small
            if (s.interrupted()) goto SUBEXPR_DOUBLE_BREAK;
            while (k != clause.end()) {
                // Count subexpressions of length 2
                std::pair<Lit, Lit> key = mkLitPair(*j, *k);

                // Add to the counter for this literal pair
                std::tr1::unordered_map<std::pair<Lit, Lit>, int>::iterator it = subexprs.find(key);
                if (it == subexprs.end()) subexprs.insert(std::make_pair(key, 1));
                else it->second++;

                k++;
            }
        }
    }
#endif
    SUBEXPR_DOUBLE_BREAK:;
    return subexprs;
}

inline std::vector< std::pair<Lit, Lit> > Solver::getFreqSubexprs(std::tr1::unordered_map<std::pair<Lit, Lit>, int>& subexpr_counts, Solver& solver, unsigned int numSubexprs) {
#define USE_QUICKSELECT_SUBEXPRS
#ifdef  USE_QUICKSELECT_SUBEXPRS
    // Copy keys
    std::vector< std::pair<Lit, Lit> > subexprs(subexpr_counts.size());
    int i = 0;
    for (std::tr1::unordered_map<std::pair<Lit, Lit>, int>::iterator it = subexpr_counts.begin(); it != subexpr_counts.end(); it++) {
        subexprs[i++] = it->first;
    }

    // Quickselect
    quickselect_count(subexprs, subexpr_counts, solver, 0, subexprs.size() - 1, numSubexprs);
    if (numSubexprs < subexprs.size()) subexprs.erase(subexprs.begin() + numSubexprs, subexprs.end());

#else
    std::vector< std::pair<Lit, Lit> > subexprs;
    std::tr1::unordered_set< std::pair<Lit, Lit> > subexprWindow;
    std::tr1::unordered_map<std::pair<Lit, Lit>, int>::iterator max = subexpr_counts.begin();

    for (unsigned int i = 0; i < numSubexprs && i < subexpr_counts.size(); i++) {
        for (std::tr1::unordered_map<std::pair<Lit, Lit>, int>::iterator it = subexpr_counts.begin(); it != subexpr_counts.end(); it++) {
            if ((subexprWindow.find(it->first) == subexprWindow.end()) && (it->second >= max->second)) {
                max = it;
            }
        }

        subexprWindow.insert(max->first);
        subexprs.push_back(max->first);
    }

#endif
    return subexprs;
}

static inline std::vector< std::tr1::unordered_set<Lit> > getLiteralSets(ClauseAllocator& ca, std::vector<CRef>& clauses) {
    // Get the set of literals for each clause
    // Time complexity: O(w k log(w k))
    std::vector< std::tr1::unordered_set<Lit> > sets;
    for (unsigned int i = 0; i < clauses.size(); i++) {
        const Clause& c = ca[clauses[i]];
        std::tr1::unordered_set<Lit> set;
        for (int j = 0; j < c.size(); j++) set.insert(c[j]);
        sets.push_back(set);
    }

    return sets;
}

std::vector< std::pair< Var, std::pair<Lit, Lit> > > Solver::user_er_add_subexpr(Solver& s, std::vector<CRef>& candidateClauses, unsigned int maxNumNewVars) {
    // Get the set of literals for each clause
    // FIXME: is there an efficient way to do this without passing in the solver's clause allocator?
    // Ideally, we want the solver object to be const when passed into this function
    std::vector< std::tr1::unordered_set<Lit> > sets = getLiteralSets(s.ca, candidateClauses);

    // Count subexpressions of length 2
    std::tr1::unordered_map<std::pair<Lit, Lit>, int> subexprs = countSubexprs(s, sets);

    // Get most frequent subexpressions
    std::vector< std::pair<Lit, Lit> > freqSubExprs = getFreqSubexprs(subexprs, s, maxNumNewVars);

    // Add extension variables
    std::vector< std::pair< Var, std::pair<Lit, Lit> > > extClauses;
    Var x = s.nVars() + s.extBuffer.size();
    for (std::vector< std::pair<Lit, Lit> >::iterator i = freqSubExprs.begin(); i != freqSubExprs.end(); i++) {
        if (!s.extVarDefs.contains(i->first, i->second)) {
            // Add extension variable
            extClauses.push_back(std::make_pair(x, *i));
            x++;
        }
    }
    return extClauses;
}

#elif ER_USER_ADD_HEURISTIC == ER_ADD_HEURISTIC_RANDOM
static inline std::vector<Var> getVarVec(ClauseAllocator& ca, std::vector<CRef>& clauses) {
    // Get set of all variables
    // Time complexity: O(w k)
    std::tr1::unordered_set<Var> vars;
    for (unsigned int i = 0; i < clauses.size(); i++)
        for (int j = 0; j < ca[clauses[i]].size(); j++)
            vars.insert(var(ca[clauses[i]][j]));
    return std::vector<Var>(vars.begin(), vars.end());
}

std::vector< std::pair< Var, std::pair<Lit, Lit> > > Solver::user_er_add_random(Solver& s, std::vector<CRef>& candidateClauses, unsigned int maxNumNewVars) {
    // Total time complexity: O(w k log(w k) + x)
    // w: window size
    // k: clause width
    // n: number of variables
    // x: maximum number of extension variables to introduce at once

    // Get set of all variables
    // Time complexity: O(w k log(w k))
    std::vector<Var> varVec = getVarVec(s.ca, candidateClauses);

    // Add extension variables
    // Time complexity: O(x)
    std::vector< std::pair< Var, std::pair<Lit, Lit> > > extClauses;
    Var x = s.nVars() + s.extBuffer.size();
    if (varVec.size() == 0) return extClauses;
    for (unsigned int i = 0; i < maxNumNewVars; i++) {
        // Sample literals at random
        int i_a = irand(s.random_seed, static_cast<int>(varVec.size()));
        int i_b = i_a;
        while (i_a == i_b) i_b = irand(s.random_seed, static_cast<int>(varVec.size()));
        Lit a = mkLit(varVec[i_a], irand(s.random_seed, 1));
        Lit b = mkLit(varVec[i_b], irand(s.random_seed, 1));

        std::pair<Lit, Lit> key = mkLitPair(a, b);
        if (!s.extVarDefs.contains(a, b)) {
            // Add extension variable
            extClauses.push_back(std::make_pair(x, key));
            x++;
        }
    }

    return extClauses;
}

#elif ER_USER_ADD_HEURISTIC == ER_ADD_HEURISTIC_GLUCOSER

std::vector< std::pair< Var, std::pair<Lit, Lit> > > Solver::user_er_add_glucosER(Solver& s, std::vector<CRef>& candidateClauses, unsigned int maxNumNewVars) {
    std::vector< std::pair< Var, std::pair<Lit, Lit> > > extDefPairs;
    // Check that we have at least 2 learnt clauses
    if (candidateClauses.size() != 2) return extDefPairs;
    Clause& c0 = s.ca[candidateClauses[0]];
    Clause& c1 = s.ca[candidateClauses[1]];

    // Check that the last 2 learnt clauses only differ by one literal
    if (c0.size() != c1.size()) return extDefPairs;
    std::tr1::unordered_set<Lit> lits;
    std::vector<Lit> diffLits;
    for (int i = 0; i < c0.size(); i++) lits.insert(c0[i]);
    for (int i = 0; i < c1.size(); i++) {
        std::tr1::unordered_set<Lit>::iterator it = lits.find(c1[i]);
        if (it == lits.end()) {
            if (diffLits.size()) return extDefPairs; // More than one different literal!
            diffLits.push_back(c1[i]);
        } else {
            lits.erase(it);
        }
    }

    // Lits and diffLits must both have size 1 here
    if (lits.size() == 1 && diffLits.size() == 1) {
        // Lits a and b must be different - otherwise, lits.size() == diffLits.size() == 0
        Lit a = *lits.begin();
        Lit b = *diffLits.begin();
        std::pair<Lit, Lit> key = mkLitPair(a, b);
        if (!s.extVarDefs.contains(a, b)) {
            // Add extension variable
            Var x = s.nVars() + s.extBuffer.size();
            extDefPairs.push_back(std::make_pair(x, key));
        }
    }

    return extDefPairs;
}

#endif

// EXTENDED RESOLUTION - variable deletion heuristic
std::tr1::unordered_set<Var> Solver::user_er_delete_all(Solver& s) {
    std::tr1::unordered_set<Var> toDelete;
    for (int i = s.originalNumVars + 1; i < s.nVars(); i++)
        toDelete.insert(i);
    return toDelete;
}

std::tr1::unordered_set<Var> Solver::user_er_delete_activity(Solver& s) {
    std::tr1::unordered_set<Var> toDelete;
    const double activityThreshold = 60; 
    for (int i = s.originalNumVars + 1; i < s.nVars(); i++) {
        if (s.activity[i] < activityThreshold) {
            toDelete.insert(i);
        }
    }
    return toDelete;
}

}