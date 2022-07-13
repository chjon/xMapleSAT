#include <algorithm>
#include <iterator>
#include <tr1/unordered_set>
#include <core/SolverER.h>
#include <core/SolverERTypes.h>

namespace Minisat {

bool SolverER::user_extFilPredicate(CRef cr) {
    const Clause& clause = solver->ca[cr];

#if ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_RANGE
    const int sz = clause.size();
    return solver->ext_min_width <= sz && sz <= solver->ext_max_width;
#elif ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LBD
    const int lbd = clause.lbd();
    return solver->ext_min_lbd <= lbd && lbd <= solver->ext_max_lbd;
#endif
}

void SolverER::user_extSelHeuristic_all(std::vector<CRef>& output, const std::vector<CRef>& input, unsigned int numClauses) {
    std::copy(input.begin(), input.end(), std::back_inserter(output));
}

void SolverER::user_extSelHeuristic_activity(std::vector<CRef>& output, const std::vector<CRef>& input, unsigned int numClauses) {
    ClauseAllocator& ca = solver->ca;
    for (const CRef cr : input) {
        CRef clauseIndex = cr, tmp = 0;
        double clauseActivity = ca[clauseIndex].activity();

        for (unsigned int i = 0; i < output.size() && i < numClauses; i++) {
            if (clauseIndex == output[i]) goto user_extSelHeuristic_activity_SKIP;
            if (clauseActivity > ca[output[i]].activity()) {
                tmp = output[i];
                output[i] = clauseIndex;
                clauseActivity = ca[tmp].activity();
                clauseIndex = tmp;
            }
        }

        if (output.size() != numClauses) output.push_back(clauseIndex);

user_extSelHeuristic_activity_SKIP: continue;
    }
}

void SolverER::user_extDefHeuristic_subexpression(std::vector<ExtDef>& extVarDefBuffer, const std::vector<CRef>& selectedClauses, unsigned int maxNumNewVars) {

}

void SolverER::user_extDefHeuristic_random(std::vector<ExtDef>& extVarDefBuffer, const std::vector<CRef>& selectedClauses, unsigned int maxNumNewVars) {
    // Total time complexity: O(w k log(w k) + x)
    // w: window size
    // k: clause width
    // n: number of variables
    // x: maximum number of extension variables to introduce at once

    // Get set of all literals
    // Time complexity: O(w k log(w k))
    ClauseAllocator& ca = solver->ca;
    std::tr1::unordered_set<Lit> litSet;
    for (unsigned int i = 0; i < selectedClauses.size(); i++) {
        Clause c = ca[selectedClauses[i]];
        for (int j = 0; j < c.size(); j++) litSet.insert(c[j]);
    }

    // NOTE:
    // size == 0 breaks randInt generation
    // size == 1 results in no valid literal pair; might as well exit here since we're checking for 0
    if (litSet.size() <= 1) return;

    std::vector<Lit> litVec = std::vector<Lit>(litSet.begin(), litSet.end());

    // Add extension variables
    // Time complexity: O(x)

    std::tr1::unordered_set< std::pair<Lit, Lit> > generatedPairs;

    // NOTE: This does a for loop instead of a while loop to avoid getting stuck
    // It's possible that we'll never be able to make as many new variables as requested if the definitions
    // already exist and there aren't enough distinct literals in the selected clauses

    Var x = solver->nVars() + extVarDefBuffer.size();
    for (unsigned int i = 0; i < maxNumNewVars; i++) {
        // Sample literals at random
        const int i_a = Solver::irand(solver->random_seed, static_cast<int>(litVec.size()));
        const int i_b = Solver::irand(solver->random_seed, static_cast<int>(litVec.size()));
        Lit a = litVec[i_a], b = litVec[i_b];
        std::pair<Lit, Lit> litPair = mkLitPair(a, b);

        // Ensure literal pair consists of different variables
        // Ensure literal pair has not already been added
        // FIXME: need to handle the case where the pair is currently queued for variable introduction in the extVarDefBuffer
        // We cannot just add directly to xdm since it might result in substitution before we introduce the new var
        if (
            var(a) == var(b) ||
            xdm.contains(a, b) ||
            generatedPairs.find(litPair) != generatedPairs.end()
        ) continue;

        // Add extension variable
        generatedPairs.insert(litPair);
        extVarDefBuffer.push_back(ExtDef { x, a, b, std::vector< std::vector<Lit> >() });
        x++;
    }
}

bool SolverER::user_extSubPredicate_size_lbd(vec<Lit>& clause) {
    if (xdm.size() == 0) return false;

#if ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_WIDTH
    const int clause_width = clause.size();
    if (clause_width < ext_sub_min_width || ext_sub_min_width > ext_sub_max_width) return false;
#endif

#if ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_LBD
    const int clause_lbd = computeLBD(clause);
    if (clause_lbd < ext_min_lbd || clause_lbd > ext_max_lbd) return false;
#endif

    return true;
}

}