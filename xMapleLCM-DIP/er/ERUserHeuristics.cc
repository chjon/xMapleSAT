/*****************************************************************************[ERUserHeuristics.cc]
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

#include <algorithm>
#include <iterator>

#include "core/ClauseDatabase.h"
#include "core/Solver.h"
#include "er/ERManager.h"
#include "er/ERTypes.h"
#include "mtl/Sort.h"

namespace Minisat {

bool ERManager::user_extFilPredicate_width(CRef cr) {
    const int sz = ca[cr].size();
    return ext_min_width <= sz && sz <= ext_max_width;
}

bool ERManager::user_extFilPredicate_lbd(CRef cr) {
    const int lbd = assignmentTrail.computeLBD(ca[cr]);
    return ext_min_lbd <= lbd && lbd <= ext_max_lbd;
}

int ERManager::numDiffs(vec<Lit>& output, const Clause& c1, const Clause& c2) {
    output.clear(); tmp_set.clear();

    // Make a set of all the literals in clause 1
    for (int i = 0; i < c1.size(); i++) tmp_set.insert(c1[i]);

    // Add all the literals in clause 2 that are not in clause 1
    for (int i = 0; i < c2.size(); i++) {
        if (tmp_set.find(c2[i]) == tmp_set.end()) {
            output.push(c2[i]);
        } else {
            tmp_set.erase(c2[i]);
        }
    }

    // Add all the literals in clause 1 that are not in clause 2
    for (Lit l : tmp_set) output.push(l);

    // Return the number of literals that are different between the two clauses
    return output.size();
}

bool ERManager::user_extFilPredicate_ler(CRef cr) {
    if (m_filteredClauses_ler.size() > 0) {
        const Clause& c1 = ca[cr];
        const Clause& c2 = ca[m_filteredClauses_ler[m_filteredClauses_ler.size() - 1]];

        if (
            c1.size() != c2.size() ||
            numDiffs(tmp_vec, c1, c2) != 2
        ) m_filteredClauses_ler.clear();
    }

    return true;
}

void ERManager::user_extSelHeuristic_all(std::vector<CRef>& output, const std::vector<CRef>& input, unsigned int numClauses) {
    const int start_index = std::max(0, static_cast<int>(input.size()) - static_cast<int>(numClauses));
    std::copy(input.begin() + start_index, input.end(), std::back_inserter(output));
}

void ERManager::user_extSelHeuristic_activity(std::vector<CRef>& output, const std::vector<CRef>& input, unsigned int numClauses) {
    for (unsigned int i = 0; i < input.size(); i++) {
        CRef clauseIndex = input[i], tmp = 0;
        double clauseActivity = ca[clauseIndex].activity();

        for (unsigned int j = 0; j < output.size() && j < numClauses; j++) {
            if (clauseIndex == output[j]) goto user_extSelHeuristic_activity_SKIP;
            if (clauseActivity > ca[output[j]].activity()) {
                tmp = output[j];
                output[j] = clauseIndex;
                clauseActivity = ca[tmp].activity();
                clauseIndex = tmp;
            }
        }

        if (output.size() != numClauses) output.push_back(clauseIndex);

user_extSelHeuristic_activity_SKIP: continue;
    }
}

// Partition elements such that clauses with larger activities are left of the pivot
// This is a helper function for quickselect
inline int partition_count(std::vector<LitPair>& db, std::tr1::unordered_map<LitPair, int>& subexpr_count, int l, int r, int pivot) {
    const int pivot_count = subexpr_count.find(db[pivot])->second;
    LitPair tmp;
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

void ERManager::quickselect_count(std::vector<LitPair>& db, std::tr1::unordered_map<LitPair, int>& subexpr_count, int l, int r, int k) {
    // Ensure we have a valid value of k
    k--;
    if (k <= 0 || k >= r - l + 1) return;
    while (!solver.interrupted()) {
        // Partition the array around last element and get position of pivot element in sorted array
        int pivot = l + randomNumberGenerator.irand(r - l + 1);
        pivot = partition_count(db, subexpr_count, l, r, pivot);

        // Update selection bounds
        if      (k < pivot) r = pivot - 1;
        else if (k > pivot) l = pivot + 1;
        else break;
    }
}

static inline std::vector<LitSet> getLiteralSets(ClauseAllocator& ca, const std::vector<CRef>& clauses) {
    // Get the set of literals for each clause
    // Time complexity: O(w k log(w k))
    std::vector<LitSet> sets;
    for (unsigned int i = 0; i < clauses.size(); i++) {
        const Clause& c = ca[clauses[i]];
        LitSet set;
        for (int j = 0; j < c.size(); j++) set.insert(c[j]);
        sets.push_back(set);
    }

    return sets;
}

inline std::vector<LitPair> ERManager::getFreqSubexprs(std::tr1::unordered_map<LitPair, int>& subexpr_counts, unsigned int numSubexprs) {
#define USE_QUICKSELECT_SUBEXPRS
#ifdef  USE_QUICKSELECT_SUBEXPRS
    // Copy keys
    std::vector<LitPair> subexprs(subexpr_counts.size());
    int i = 0;
    for (std::tr1::unordered_map<LitPair, int>::iterator it = subexpr_counts.begin(); it != subexpr_counts.end(); it++) {
        subexprs[i++] = it->first;
    }

    // Quickselect
    quickselect_count(subexprs, subexpr_counts, 0, subexprs.size() - 1, numSubexprs);
    if (numSubexprs < subexprs.size()) subexprs.erase(subexprs.begin() + numSubexprs, subexprs.end());

#else
    std::vector<LitPair> subexprs;
    std::tr1::unordered_set<LitPair> subexprWindow;
    std::tr1::unordered_map<LitPair, int>::iterator max = subexpr_counts.begin();

    for (unsigned int i = 0; i < numSubexprs && i < subexpr_counts.size(); i++) {
        for (std::tr1::unordered_map<LitPair, int>::iterator it = subexpr_counts.begin(); it != subexpr_counts.end(); it++) {
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

static inline void addIntersectionToSubexprs(std::tr1::unordered_map<LitPair, int>& subexprs, const std::vector<Lit>& intersection) {
    // Time complexity: O(k^2)
    for (unsigned int i = 0; i < intersection.size(); i++) {
        for (unsigned int j = i + 1; j < intersection.size(); j++) {
            // Count subexpressions of length 2
            LitPair key = mkLitPair(intersection[i], intersection[j]);

            // Add to the counter for this literal pair
            std::tr1::unordered_map<LitPair, int>::iterator it = subexprs.find(key);
            if (it == subexprs.end()) subexprs.insert(std::make_pair(key, 1));
            else it->second++;
        }
    }
}

inline std::tr1::unordered_map<LitPair, int> ERManager::countSubexprs(std::vector<LitSet>& sets) {
    // Count subexpressions by looking at intersections
    // Time complexity: O(w k^2)
    std::tr1::unordered_map<LitPair, int> subexprs;
#if ER_USER_ADD_SUBEXPR_SET_INTERSECTION
    for (unsigned int i = 0; i < sets.size(); i++) {
        for (unsigned int j = i + 1; j < sets.size(); j++) {
            // TODO: Check if we've already processed a pair of clauses (cache) and add their counts

            // We might spend a lot of time here - exit if interrupted
            // FIXME: ideally, we wouldn't have to check for this at all if the sets of literals were sufficiently small
            if (solver.interrupted()) goto SUBEXPR_DOUBLE_BREAK;
            std::vector<Lit> intersection(sets[i].size() + sets[j].size());
            std::vector<Lit>::iterator it = std::set_intersection(sets[i].begin(), sets[i].end(), sets[j].begin(), sets[j].end(), intersection.begin());
            intersection.resize(it - intersection.begin());
            addIntersectionToSubexprs(subexprs, intersection);
        }
    }
#else
    // Naive algo: for each clause (m), increment a counter for every pair of lits in the clause (k^2), then sort the pairs (n^2)
    // O(m * k^2 + n^2)
    // this is very bad if the clauses are large

    for (unsigned int i = 0; i < sets.size(); i++) {
        LitSet& clause = sets[i];
        for (LitSet::iterator j = clause.begin(); j != clause.end(); j++) {
            LitSet::iterator k = j; k++;
            // We might spend a lot of time here - exit if interrupted
            // FIXME: ideally, we wouldn't have to check for this at all if the sets of literals were sufficiently small
            if (solver.interrupted()) goto SUBEXPR_DOUBLE_BREAK;
            while (k != clause.end()) {
                // Count subexpressions of length 2
                LitPair key = mkLitPair(*j, *k);

                // Add to the counter for this literal pair
                std::tr1::unordered_map<LitPair, int>::iterator it = subexprs.find(key);
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

void ERManager::user_extDefHeuristic_subexpression(std::vector<ExtDef>& extVarDefBuffer, const std::vector<CRef>& selectedClauses, unsigned int maxNumNewVars) {
    // Get the set of literals for each clause
    // FIXME: is there an efficient way to do this without passing in the solver's clause allocator?
    // Ideally, we want the solver object to be const when passed into this function
    std::vector<LitSet> sets = getLiteralSets(ca, selectedClauses);

    // Count subexpressions of length 2
    std::tr1::unordered_map<LitPair, int> subexprs = countSubexprs(sets);

    // Get most frequent subexpressions
    std::vector<LitPair> freqSubExprs = getFreqSubexprs(subexprs, maxNumNewVars);

    // Add extension variables
    std::tr1::unordered_set<LitPair> generatedPairs;
    Var x = assignmentTrail.nVars() + extVarDefBuffer.size();
    for (std::vector<LitPair>::iterator i = freqSubExprs.begin(); i != freqSubExprs.end(); i++) {
        Lit a = i->first, b = i->second;
        if (isValidDefPair(a, b, generatedPairs)) {
            // Add extension variable
            generatedPairs.insert(*i);
            extVarDefBuffer.push_back(ExtDef { mkLit(x), a, b, std::vector< std::vector<Lit> >() });
            x++;
        }
    }
}

void ERManager::user_extDefHeuristic_random(std::vector<ExtDef>& extVarDefBuffer, const std::vector<CRef>& selectedClauses, unsigned int maxNumNewVars) {
    // Total time complexity: O(w k log(w k) + x)
    // w: window size
    // k: clause width
    // n: number of variables
    // x: maximum number of extension variables to introduce at once

    // Get set of all literals
    // Time complexity: O(w k log(w k))
    LitSet litSet;
    for (unsigned int i = 0; i < selectedClauses.size(); i++) {
        Clause& c = ca[selectedClauses[i]];
        for (int j = 0; j < c.size(); j++) litSet.insert(c[j]);
    }

    // NOTE:
    // size == 0 breaks randInt generation
    // size == 1 results in no valid literal pair; might as well exit here since we're checking for 0
    if (litSet.size() <= 1) return;

    std::vector<Lit> litVec = std::vector<Lit>(litSet.begin(), litSet.end());

    // Add extension variables
    // Time complexity: O(x)

    std::tr1::unordered_set<LitPair> generatedPairs;

    // NOTE: This does a for loop instead of a while loop to avoid getting stuck
    // It's possible that we'll never be able to make as many new variables as requested if the definitions
    // already exist and there aren't enough distinct literals in the selected clauses

    const unsigned int MAX_RETRIES = 10;
    Var x = assignmentTrail.nVars() + extVarDefBuffer.size();
    for (unsigned int i = 0; i < maxNumNewVars; i++) {
        // Sample literals at random
        for (unsigned int j = 0; j < MAX_RETRIES; j++) {
            const int i_a = randomNumberGenerator.irand(static_cast<int>(litVec.size()));
            const int i_b = randomNumberGenerator.irand(static_cast<int>(litVec.size()));
            Lit a = litVec[i_a], b = litVec[i_b];
            LitPair litPair = mkLitPair(a, b);

            if (isValidDefPair(a, b, generatedPairs)) {
                // Add extension variable
                generatedPairs.insert(litPair);
                extVarDefBuffer.push_back(ExtDef { mkLit(x), a, b, std::vector< std::vector<Lit> >() });
                x++;
                break;
            }
        }
    }
}

void ERManager::user_extDefHeuristic_ler(std::vector<ExtDef>& extVarDefBuffer, const std::vector<CRef>& selectedClauses, unsigned int maxNumNewVars) {
    // Add extension variables
    std::tr1::unordered_set<LitPair> generatedPairs;
    Var x = assignmentTrail.nVars() + extVarDefBuffer.size();
    for (unsigned int i = 1; i < selectedClauses.size(); i++) {
        int n = numDiffs(tmp_vec, ca[selectedClauses[i]], ca[selectedClauses[i - 1]]);
        assert(n == 2);
        const Lit a = ~tmp_vec[0], b = ~tmp_vec[1];
        LitPair litPair = mkLitPair(a, b);
        if (isValidDefPair(a, b, generatedPairs)) {
            // Add extension variable
            generatedPairs.insert(litPair);
            extVarDefBuffer.push_back(ExtDef { mkLit(x), a, b, std::vector< std::vector<Lit> >() });
            x++;
        }
    }
}

bool ERManager::user_extSubPredicate_size_lbd(vec<Lit>& clause) {
    if (xdm.size() == 0) return false;

#if ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_WIDTH
    const int clause_width = clause.size();
    if (clause_width < ext_sub_min_width) return false;
#endif

#if ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_LBD
    const int clause_lbd = assignmentTrail.computeLBD(clause);
    if (clause_lbd < ext_min_lbd || clause_lbd > ext_max_lbd) return false;
#endif

    return true;
}

void ERManager::user_extDelPredicateSetup_none() {}
void ERManager::user_extDelPredicateSetup_activity() {
#if ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY
    m_threshold_activity = ext_act_threshold;
#elif ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY2
    // Copy activities for current variables
    const vec<double>& activity = solver.branchingHeuristicManager.getActivity();
    vec<double> currentActivity(originalNumVars);
    for (int i = 0; i < originalNumVars; i++) currentActivity.push(activity[i]);
    for (auto it = extDefs.begin(); it != extDefs.end(); it++) currentActivity.push(activity[it->first]);
    
    // Compute threshold activity
    sort(currentActivity);
    m_threshold_activity = currentActivity[(int)((currentActivity.size() - 1) * ext_act_threshold)];
#endif
}

bool ERManager::user_extDelPredicate_none(Var x) { return false; }
bool ERManager::user_extDelPredicate_all(Var x) { return true; }
bool ERManager::user_extDelPredicate_activity(Var x) {
#if ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY || ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY2
    const double activity = solver.branchingHeuristicManager.getActivity()[x];
    return activity < m_threshold_activity;
#else
    return false;
#endif
}

}
