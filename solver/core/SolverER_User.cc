#include <algorithm>
#include "core/Solver.h"

namespace Minisat {

static inline std::pair<Lit, Lit> mkLitPair(Lit a, Lit b) {
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

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
    // FIXME: this is probably inefficient, but there isn't a preexisting data structure which keeps these in sorted order
    // Optimization idea: sort each of the clause DBs and then pick the top k (could also use heap sort)

    // Optimization idea: use a quicksort-type of algo for expected-linear time
    std::vector<CRef> clauseWindow;
    for (int i = 0; i < s.nClauses   (); i++) addClauseToWindow(s.ca, clauseWindow, s.clauses   [i], numClauses);
    for (int i = 0; i < s.nLearnts   (); i++) addClauseToWindow(s.ca, clauseWindow, s.learnts   [i], numClauses);
    for (int i = 0; i < s.nExtLearnts(); i++) addClauseToWindow(s.ca, clauseWindow, s.extLearnts[i], numClauses);
    for (int i = 0; i < s.nExtDefs   (); i++) addClauseToWindow(s.ca, clauseWindow, s.extDefs   [i], numClauses);
    return clauseWindow;
}

// Partition elements such that clauses with smaller activities are left of the pivot
// This is a helper function for quickselect
static inline int partition(vec<CRef>& db, ClauseAllocator& ca, int l, int r, int pivot) {
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
void Solver::quickselect(vec<CRef>& db, Solver& solver, int l, int r, int k) {
    // Ensure we have a valid value of k
    if (k <= 0 || k >= r - l + 1) return;
    while (!solver.interrupted()) {
        // Partition the array around last element and get position of pivot element in sorted array
        int pivot = l + solver.irand(solver.random_seed, r - l + 1);
        pivot = partition(db, solver.ca, l, r, pivot);

        if (k < pivot) {
            r = pivot - 1;
        } else if (k > pivot) {
            l = pivot + 1;
        } else {
            break;
        }
    }
}

// Get the elements with the k largest activities and store them in the target vector
void Solver::copyKLargest(vec<CRef>& target, vec<CRef>& source, Solver& solver, unsigned int k) {
#define ER_ALLOW_DB_MODIFICATION false
#if ER_ALLOW_DB_MODIFICATION
    quickselect(source, solver, 0, source.size() - 1, k);
    for (unsigned int i = 0; i < std::min(static_cast<unsigned int>(source.size()), k); i++) target.push(source[i]);
#else
    // TODO: is making this copy necessary?
    vec<CRef> copy;
    source.copyTo(copy);
    quickselect(copy, solver, 0, copy.size() - 1, k);
    for (unsigned int i = 0; i < std::min(static_cast<unsigned int>(copy.size()), k); i++) target.push(copy[i]);
#endif
}

std::vector<CRef> Solver::user_er_select_activity2(Solver& s, unsigned int numClauses) {
    // Find the variables in the clauses with the top k highest activities
    // Uses a quicksort-type of algo for expected-linear time
    // Adapted from: https://en.wikipedia.org/wiki/Quickselect

    vec<CRef> clauses;
    copyKLargest(clauses, s.clauses   , s, numClauses);
    copyKLargest(clauses, s.learnts   , s, numClauses);
    copyKLargest(clauses, s.extLearnts, s, numClauses);
    copyKLargest(clauses, s.extDefs   , s, numClauses);

    std::vector<CRef> clauseWindow;
    quickselect(clauses, s, 0, clauses.size() - 1, numClauses);
    for (unsigned int i = 0; i < std::min(static_cast<unsigned int>(clauses.size()), numClauses); i++) clauseWindow.push_back(clauses[i]);
    return clauseWindow;
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

static inline std::tr1::unordered_map<std::pair<Lit, Lit>, int> countSubexprs(const Solver& s, std::vector< std::tr1::unordered_set<Lit> >& sets) {
    // Count subexpressions by looking at intersections
    // Time complexity: O(w k^2)
    std::tr1::unordered_map<std::pair<Lit, Lit>, int> subexprs;

// #define USE_SET_INTERSECTION
#ifdef USE_SET_INTERSECTION
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
    SUBEXPR_DOUBLE_BREAK:;
#else
    for (unsigned int i = 0; i < sets.size() && !s.interrupted(); i++) {
        std::tr1::unordered_set<Lit>& clause = sets[i];
        // printf("%u/%lu: %lu\n", i, sets.size(), clause.size());
        if (clause.size() > 100) continue;

        for (std::tr1::unordered_set<Lit>::iterator j = clause.begin(); j != clause.end(); j++) {
            std::tr1::unordered_set<Lit>::iterator k = j; k++;
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
    return subexprs;
}

static inline std::tr1::unordered_set< std::pair<Lit, Lit> > getFreqSubexprs(std::tr1::unordered_map<std::pair<Lit, Lit>, int>& subexprs, unsigned int numSubexprs) {
    std::tr1::unordered_set< std::pair<Lit, Lit> > subexprWindow;
    std::tr1::unordered_map<std::pair<Lit, Lit>, int>::iterator max = subexprs.begin();
    for (unsigned int i = 0; i < numSubexprs && i < subexprs.size(); i++) {
        for (std::tr1::unordered_map<std::pair<Lit, Lit>, int>::iterator it = subexprs.begin(); it != subexprs.end(); it++) {
            if ((subexprWindow.find(it->first) == subexprWindow.end()) && (it->second >= max->second)) {
                max = it;
            }
        }

        subexprWindow.insert(max->first);
    }

    return subexprWindow;
}

static inline bool contains (std::tr1::unordered_map< Lit, std::tr1::unordered_map<Lit, Lit> > extVarMap, Lit a, Lit b) {
    std::tr1::unordered_map< Lit, std::tr1::unordered_map<Lit, Lit> >::const_iterator it1 = extVarMap.find(a);
    if (it1 == extVarMap.end()) return false;
    std::tr1::unordered_map<Lit, Lit>::const_iterator it2 = it1->second.find(b);
    return it2 != it1->second.end();
}

// EXTENDED RESOLUTION - variable definition heuristic
std::tr1::unordered_map< Var, std::pair<Lit, Lit> > Solver::user_er_add_subexpr(Solver& s, std::vector<CRef>& candidateClauses, unsigned int maxNumNewVars) {
    // Get the set of literals for each clause
    // FIXME: is there an efficient way to do this without passing in the solver's clause allocator?
    // Ideally, we want the solver object to be const when passed into this function
    std::vector< std::tr1::unordered_set<Lit> > sets = getLiteralSets(s.ca, candidateClauses);

    // Count subexpressions of length 2
    std::tr1::unordered_map<std::pair<Lit, Lit>, int> subexprs = countSubexprs(s, sets);

    // Get most frequent subexpressions
    std::tr1::unordered_set< std::pair<Lit, Lit> > freqSubExprs = getFreqSubexprs(subexprs, maxNumNewVars);

    // Add extension variables
    std::tr1::unordered_map< Var, std::pair<Lit, Lit> > extClauses;
    Var x = s.nVars();
    for (std::tr1::unordered_set< std::pair<Lit, Lit> >::iterator i = freqSubExprs.begin(); i != freqSubExprs.end(); i++) {
        if (!contains(s.extVarDefs, i->first, i->second)) {
            // Add extension variable
            extClauses.insert(std::make_pair(x, *i));
            x++;
        }
    }
    return extClauses;
}

static inline std::vector<Var> getVarVec(ClauseAllocator& ca, std::vector<CRef>& clauses) {
    // Get set of all variables
    // Time complexity: O(w k log(w k))
    std::tr1::unordered_set<Var> vars;
    for (unsigned int i = 0; i < clauses.size(); i++)
        for (int j = 0; j < ca[clauses[i]].size(); j++)
            vars.insert(var(ca[clauses[i]][j]));
    return std::vector<Var>(vars.begin(), vars.end());
}

// EXTENDED RESOLUTION - variable definition heuristic
std::tr1::unordered_map< Var, std::pair<Lit, Lit> > Solver::user_er_add_random(Solver& s, std::vector<CRef>& candidateClauses, unsigned int maxNumNewVars) {
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
    std::tr1::unordered_map< Var, std::pair<Lit, Lit> > extClauses;
    Var x = s.nVars();
    for (unsigned int i = 0; i < maxNumNewVars; i++) {
        // Sample literals at random
        int i_a = irand(s.random_seed, static_cast<int>(varVec.size()));
        int i_b = i_a;
        while (i_a == i_b) i_b = irand(s.random_seed, static_cast<int>(varVec.size()));
        Lit a = mkLit(varVec[i_a], irand(s.random_seed, 1));
        Lit b = mkLit(varVec[i_b], irand(s.random_seed, 1));

        std::pair<Lit, Lit> key = mkLitPair(a, b);
        if (!contains(s.extVarDefs, a, b)) {
            // Add extension variable
            extClauses.insert(std::make_pair(x, key));
            x++;
        }
    }

    return extClauses;
}

// EXTENDED RESOLUTION - variable deletion heuristic
std::vector<Var> Solver::user_er_delete_all(Solver& s) {
    std::vector<Var> toDelete;
    for (int i = s.originalNumVars + 1; i < s.nVars(); i++)
        toDelete.push_back(i);
    return toDelete;
}

}