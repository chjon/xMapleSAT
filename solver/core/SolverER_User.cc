#include <algorithm>
#include <iterator>
#include <core/SolverER.h>
#include <core/SolverERTypes.h>

namespace Minisat {

bool SolverER::user_extFilPredicate_width(CRef cr) {
    const int sz = solver->ca[cr].size();
    return ext_min_width <= sz && sz <= ext_max_width;
}

bool SolverER::user_extFilPredicate_lbd(CRef cr) {
    const int lbd = solver->lbd(solver->ca[cr]);
    return ext_min_lbd <= lbd && lbd <= ext_max_lbd;
}

void SolverER::user_extSelHeuristic_all(std::vector<CRef>& output, const std::vector<CRef>& input, unsigned int numClauses) {
    std::copy(input.begin(), input.end(), std::back_inserter(output));
}

void SolverER::user_extSelHeuristic_activity(std::vector<CRef>& output, const std::vector<CRef>& input, unsigned int numClauses) {
    ClauseAllocator& ca = solver->ca;
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

void SolverER::quickselect_count(std::vector<LitPair>& db, std::tr1::unordered_map<LitPair, int>& subexpr_count, int l, int r, int k) {
    // Ensure we have a valid value of k
    k--;
    if (k <= 0 || k >= r - l + 1) return;
    while (!solver->asynch_interrupt) {
        // Partition the array around last element and get position of pivot element in sorted array
        int pivot = l + solver->irand(solver->random_seed, r - l + 1);
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

inline std::vector<LitPair> SolverER::getFreqSubexprs(std::tr1::unordered_map<LitPair, int>& subexpr_counts, unsigned int numSubexprs) {
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

inline std::tr1::unordered_map<LitPair, int> SolverER::countSubexprs(std::vector<LitSet>& sets) {
    // Count subexpressions by looking at intersections
    // Time complexity: O(w k^2)
    std::tr1::unordered_map<LitPair, int> subexprs;
#if ER_USER_ADD_SUBEXPR_SET_INTERSECTION
    for (unsigned int i = 0; i < sets.size(); i++) {
        for (unsigned int j = i + 1; j < sets.size(); j++) {
            // TODO: Check if we've already processed a pair of clauses (cache) and add their counts

            // We might spend a lot of time here - exit if interrupted
            // FIXME: ideally, we wouldn't have to check for this at all if the sets of literals were sufficiently small
            if (solver->asynch_interrupt) goto SUBEXPR_DOUBLE_BREAK;
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
            if (solver->asynch_interrupt) goto SUBEXPR_DOUBLE_BREAK;
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

void SolverER::user_extDefHeuristic_subexpression(std::vector<ExtDef>& extVarDefBuffer, const std::vector<CRef>& selectedClauses, unsigned int maxNumNewVars) {
    // Get the set of literals for each clause
    // FIXME: is there an efficient way to do this without passing in the solver's clause allocator?
    // Ideally, we want the solver object to be const when passed into this function
    std::vector<LitSet> sets = getLiteralSets(solver->ca, selectedClauses);

    // Count subexpressions of length 2
    std::tr1::unordered_map<LitPair, int> subexprs = countSubexprs(sets);

    // Get most frequent subexpressions
    std::vector<LitPair> freqSubExprs = getFreqSubexprs(subexprs, maxNumNewVars);

    // Add extension variables
    std::tr1::unordered_set<LitPair> generatedPairs;
    Var x = solver->nVars() + extVarDefBuffer.size();
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

void SolverER::user_extDefHeuristic_random(std::vector<ExtDef>& extVarDefBuffer, const std::vector<CRef>& selectedClauses, unsigned int maxNumNewVars) {
    // Total time complexity: O(w k log(w k) + x)
    // w: window size
    // k: clause width
    // n: number of variables
    // x: maximum number of extension variables to introduce at once

    // Get set of all literals
    // Time complexity: O(w k log(w k))
    LitSet litSet;
    for (unsigned int i = 0; i < selectedClauses.size(); i++) {
        Clause& c = solver->ca[selectedClauses[i]];
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
    Var x = solver->nVars() + extVarDefBuffer.size();
    for (unsigned int i = 0; i < maxNumNewVars; i++) {
        // Sample literals at random
        for (unsigned int j = 0; j < MAX_RETRIES; j++) {
            const int i_a = Solver::irand(solver->random_seed, static_cast<int>(litVec.size()));
            const int i_b = Solver::irand(solver->random_seed, static_cast<int>(litVec.size()));
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

bool SolverER::user_extSubPredicate_size_lbd(vec<Lit>& clause) {
    if (xdm.size() == 0) return false;

#if ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_WIDTH
    const int clause_width = clause.size();
    if (clause_width < ext_sub_min_width) return false;
#endif

#if ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_LBD
    const int clause_lbd = clause.lbd();
    if (clause_lbd < ext_min_lbd || clause_lbd > ext_max_lbd) return false;
#endif

    return true;
}

bool SolverER::user_extDelPredicate_none(Var x) { return false; }
bool SolverER::user_extDelPredicate_all(Var x) { return true; }
bool SolverER::user_extDelPredicate_activity(Var x) {
    return solver->activity[x] < ext_act_threshold;
}

}