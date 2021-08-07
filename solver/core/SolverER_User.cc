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

static inline void addIntersectionToSubexprs(std::map<std::set<Lit>, int>& subexprs, const std::vector<Lit>& intersection) {
    std::set<Lit> subexpr;
    for (unsigned int i = 0; i < intersection.size(); i++) {
        for (unsigned int j = i + 1; j < intersection.size(); j++) {
            // Count subexpressions of length 2
            subexpr.clear();
            subexpr.insert(intersection[i]);
            subexpr.insert(intersection[j]);

            // Add to the counter for this literal pair
            std::map<std::set<Lit>, int>::iterator it = subexprs.find(subexpr);
            if (it == subexprs.end()) subexprs.insert(std::make_pair(subexpr, 1));
            else it->second++;
        }
    }
}

static inline std::vector< std::set<Lit> > getLiteralSets(ClauseAllocator& ca, std::vector<CRef>& clauses) {
    // Get the set of literals for each clause
    std::vector< std::set<Lit> > sets;
    for (unsigned int i = 0; i < clauses.size(); i++) {
        const Clause& c = ca[clauses[i]];
        std::set<Lit> set;
        for (int j = 0; j < c.size(); j++) set.insert(c[j]);
        sets.push_back(set);
    }

    return sets;
}

static inline std::map<std::set<Lit>, int> countSubexprs(const Solver& s, std::vector< std::set<Lit> >& sets) {
    // Count subexpressions by looking at intersections
    std::map<std::set<Lit>, int> subexprs;
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
    SUBEXPR_DOUBLE_BREAK:
    return subexprs;
}

static inline std::set< std::set<Lit> > getFreqSubexprs(std::map<std::set<Lit>, int>& subexprs, unsigned int numSubexprs) {
    std::set< std::set<Lit> > subexprWindow;
    std::map<std::set<Lit>, int>::iterator max = subexprs.begin();
    for (unsigned int i = 0; i < numSubexprs && i < subexprs.size(); i++) {
        for (std::map<std::set<Lit>, int>::iterator it = subexprs.begin(); it != subexprs.end(); it++) {
            if ((subexprWindow.find(it->first) == subexprWindow.end()) && (it->second >= max->second)) {
                max = it;
            }
        }

        subexprWindow.insert(max->first);
    }

    return subexprWindow;
}

// EXTENDED RESOLUTION - variable definition heuristic
std::map< Var, std::pair<Lit, Lit> > Solver::user_er_add_subexpr(Solver& s, std::vector<CRef>& candidateClauses, unsigned int maxNumNewVars) {
    // Get the set of literals for each clause
    // FIXME: is there an efficient way to do this without passing in the solver's clause allocator?
    // Ideally, we want the solver object to be const when passed into this function
    std::vector< std::set<Lit> > sets = getLiteralSets(s.ca, candidateClauses);

    // Count subexpressions of length 2
    std::map<std::set<Lit>, int> subexprs = countSubexprs(s, sets);

    // Get most frequent subexpressions
    std::set< std::set<Lit> > freqSubExprs = getFreqSubexprs(subexprs, maxNumNewVars);

    // Add extension variables
    std::map< Var, std::pair<Lit, Lit> > extClauses;
    Var x = s.nVars();
    for (std::set< std::set<Lit> >::iterator i = freqSubExprs.begin(); i != freqSubExprs.end(); i++) {
        std::set<Lit>::const_iterator it = i->begin();
        Lit a = *it; it++;
        Lit b = *it; it++;

        std::pair<Lit, Lit> key = mkLitPair(a, b);
        if (s.extVarDefs.find(key) == s.extVarDefs.end()) {
            // Add extension variable
            extClauses.insert(std::make_pair(x, key));
            x++;
        }
    }
    return extClauses;
}

static inline std::vector<Var> getVarVec(ClauseAllocator& ca, std::vector<CRef>& clauses) {
    // Get set of all variables
    std::set<Var> vars;
    for (unsigned int i = 0; i < clauses.size(); i++)
        for (int j = 0; j < ca[clauses[i]].size(); j++)
            vars.insert(var(ca[clauses[i]][j]));
    return std::vector<Var>(vars.begin(), vars.end());
}

// EXTENDED RESOLUTION - variable definition heuristic
std::map< Var, std::pair<Lit, Lit> > Solver::user_er_add_random(Solver& s, std::vector<CRef>& candidateClauses, unsigned int maxNumNewVars) {
    // Get set of all variables
    std::vector<Var> varVec = getVarVec(s.ca, candidateClauses);

    // Add extension variables
    std::map< Var, std::pair<Lit, Lit> > extClauses;
    Var x = s.nVars();
    for (unsigned int i = 0; i < maxNumNewVars; i++) {
        // Sample literals at random
        int i_a = irand(s.random_seed, static_cast<int>(varVec.size()));
        int i_b = i_a;
        while (i_a == i_b) i_b = irand(s.random_seed, static_cast<int>(varVec.size()));
        Lit a = mkLit(varVec[i_a], irand(s.random_seed, 1));
        Lit b = mkLit(varVec[i_b], irand(s.random_seed, 1));

        std::pair<Lit, Lit> key = mkLitPair(a, b);
        if (s.extVarDefs.find(key) == s.extVarDefs.end()) {
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