#include "core/Solver.h"
#include "mtl/Sort.h"

// Template specializations for hashing
namespace std { namespace tr1 {
    template<>
    std::size_t std::tr1::hash<std::pair<Minisat::Lit, Minisat::Lit> >::operator()(std::pair<Minisat::Lit, Minisat::Lit> p) const {
        return std::size_t(p.first.x) << 32 | p.second.x;
    }

    template<>
    std::size_t std::tr1::hash<Minisat::Lit>::operator()(Minisat::Lit p) const {
        return p.x;
    }
}}

namespace Minisat {

static inline std::pair<Lit, Lit> mkLitPair(Lit a, Lit b) {
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

static inline void removeLits(std::tr1::unordered_set<Lit>& set, const std::vector<Lit>& toRemove) {
    for (std::vector<Lit>::const_iterator it = toRemove.begin(); it != toRemove.end(); it++) set.erase(*it);
}

inline void Solver::er_substitute(vec<Lit>& clause, struct LitPairMap& extVarDefs) {
    // Get set of all literals in clause
    std::tr1::unordered_set<Lit> learntLits;
    for (int i = 1; i < clause.size(); i++) learntLits.insert(clause[i]);

    // Possible future investigation: should we ignore really long clauses in order to save time?
    // This would need a command-line option to toggle ignoring the clauses

    // Check for extension variables over literal pairs
    for (int i = 1; i < clause.size(); i++) {
        // Check if any extension variables are defined over this literal
        std::tr1::unordered_map< Lit, std::tr1::unordered_map<Lit, Lit> >::const_iterator tmp1 = extVarDefs.map.find(clause[i]);
        if (tmp1 == extVarDefs.map.end()) continue;

        // TODO: would it be better to iterate over the extension definitions here?
        // Also consider sorting or doing a set intersection here
        const std::tr1::unordered_map<Lit, Lit>& possibleDefs = tmp1->second;
        for (int j = i + 1; j < clause.size(); j++) {
            // Check if any extension variables are defined over this literal pair
            std::tr1::unordered_map<Lit, Lit>::const_iterator tmp2 = possibleDefs.find(clause[j]);
            if (tmp2 == possibleDefs.end()) continue;

            // Replace disjunction with intersection
            learntLits.erase(tmp1->first);
            learntLits.erase(tmp2->first);
            learntLits.insert(tmp2->second);
            goto ER_SUBSTITUTE_DOUBLE_BREAK;
        }
    }
    ER_SUBSTITUTE_DOUBLE_BREAK:;

    // Generate reduced learnt clause
    std::tr1::unordered_set<Lit>::iterator it = learntLits.begin();
    unsigned int i;
    for (i = 1; i <= learntLits.size(); i++) {
        clause[i] = *it;
        it++;
    }
    clause.shrink(clause.size() - i);
}

void Solver::substituteExt(vec<Lit>& out_learnt) {
    extTimerStart();
    // Ensure we have extension variables
    if (nExtVars() > 0) {
#if ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_WIDTH
        // Check clause width
        int clause_width = out_learnt.size();
        if (clause_width >= ext_sub_min_width &&
            clause_width <= ext_sub_max_width
        )
#endif
        {
#if ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_LBD
            // Check LBD
            int clause_lbd = lbd(out_learnt);
            if (clause_lbd >= ext_min_lbd && clause_lbd <= ext_max_lbd)
#endif
            {
                er_substitute(out_learnt, extVarDefs);
            }
        }
    }

    extTimerStop(ext_sub_overhead);
}

// Prioritize branching on a given set of literals
inline void Solver::er_prioritize(const std::vector<Var>& toPrioritize) {
    const double desiredActivity = activity[order_heap[0]] * 1.5;
    for (std::vector<Var>::const_iterator i = toPrioritize.begin(); i != toPrioritize.end(); i++) {
        Var v = *i;
        // Prioritize branching on our extension variables
        activity[v] = desiredActivity;
#if EXTENSION_FORCE_BRANCHING
        // This forces branching because of how branching is implemented when ANTI_EXPLORATION is turned on
        // FIXME: this only forces branching on the last extension variable we add here - maybe add a queue for force branch variables?
        canceled[v] = conflicts;
#endif
        if (order_heap.inHeap(v)) order_heap.decrease(v);
    }
}

inline std::vector<Var> Solver::er_add(
    vec<CRef>& er_def_db,
    struct LitPairMap& er_def_map,
    const std::vector< std::pair< Var, std::pair<Lit, Lit> > >& newDefMap
) {
    // Add extension variables
    // It is the responsibility of the user heuristic to ensure that we do not have pre-existing extension variables
    // for the provided literal pairs

    // TODO: don't add the extension variable in the case where we have x1 = (a v b) and x2 = (x1 v a)
    // TODO: don't add the extension variable in the case where we have x1 = (a v b) and x2 = (x1 v -a)
    std::vector<Var> new_variables;
    for (unsigned int i = 0; i < newDefMap.size(); i++) {
        new_variables.push_back(newVar());
    }

    // Add extension clauses
    for (std::vector< std::pair< Var, std::pair<Lit, Lit> > >::const_iterator i = newDefMap.begin(); i != newDefMap.end(); i++) {
        // Get literals
        Lit x = mkLit(i->first);
        Lit a = i->second.first;
        Lit b = i->second.second;
        assert(var(x) > var(a) && var(x) > var(b));

        // Create extension clauses and add them to the extension definition database
        addClauseToDB(er_def_db, ~x, a, b);
        addClauseToDB(er_def_db, x, ~a);
        addClauseToDB(er_def_db, x, ~b);

        // Save definition
        er_def_map.insert(x, a, b);
    }

    return new_variables;
}

void Solver::generateExtVars (
    std::vector<CRef>(*er_select_heuristic)(Solver&, unsigned int),
    std::vector< std::pair< Var, std::pair<Lit, Lit> > >(*er_add_heuristic)(Solver&, std::vector<CRef>&, unsigned int),
    unsigned int numClausesToConsider,
    unsigned int maxNumNewVars
) {
    // Get extension clauses according to heuristic
    extTimerStart();
    std::vector<CRef> candidateClauses = er_select_heuristic(*this, numClausesToConsider);
    extTimerStop(ext_sel_overhead);

    // Get extension variables according to heuristic
    extTimerStart();
    const std::vector< std::pair< Var, std::pair<Lit, Lit> > > newDefMap = er_add_heuristic(*this, candidateClauses, maxNumNewVars);
    extBuffer.insert(extBuffer.end(), newDefMap.begin(), newDefMap.end());
    extTimerStop(ext_add_overhead);
}

// Add extension variables to our data structures and prioritize branching on them.
// This calls a heuristic function which is responsible for identifying extension variable
// definitions and adding the appropriate clauses and variables.
void Solver::addExtVars() {
    // Add the extension variables to our data structures
    extTimerStart();
    const std::vector<Var> new_variables = er_add(extDefs, extVarDefs, extBuffer);
    extBuffer.clear();
    er_prioritize(new_variables);
    extTimerStop(ext_add_overhead);
}

void Solver::delExtVars(Minisat::vec<Minisat::CRef>& db, const std::tr1::unordered_set<Var>& varsToDeleteSet) {
    int i, j;

    // Delete clauses which contain the extension variable
    // TODO: is there a more efficient way to implement this? e.g. have a list of clauses for each extension variable?
    for (i = j = 0; i < db.size(); i++) {
        Clause& c = ca[db[i]];
        bool containsVarToDelete = false;
        for (int k = 0; k < c.size(); k++) {
            if (varsToDeleteSet.find(var(c[k])) != varsToDeleteSet.end()) {
                containsVarToDelete = true;
                break;
            }
        }

        if (!locked(c) && containsVarToDelete) {
#if ER_USER_FILTER_HEURISTIC != ER_FILTER_HEURISTIC_NONE
            user_er_filter_delete_incremental(db[i]);
#endif
            removeClause(db[i]);
        } else {
            db[j++] = db[i];
        }
    }
    db.shrink(i - j);

#if ER_USER_FILTER_HEURISTIC != ER_FILTER_HEURISTIC_NONE
    user_er_filter_delete_flush();
#endif
}

void Solver::delExtVars(std::tr1::unordered_set<Var>(*er_delete_heuristic)(Solver&)) {
    extTimerStart();

    // Delete variables

    // TODO: add an option to switch between these two modes

#if ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ALL
    // option 1: delete all clauses containing the extension variables
    extLearnts.clear();
    extDefs.clear(); 
    extVarDefs.clear();
#elif ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY
    // Get variables to delete
    std::tr1::unordered_set<Var> varsToDelete = er_delete_heuristic(*this);
    delExtVars(extLearnts, varsToDelete);
    delExtVars(extDefs, varsToDelete);   

    // Remove definitions from data structures
    extVarDefs.erase(varsToDelete);
#endif
    // option 2: substitute extension variable with definition
    // TODO

    extTimerStop(ext_delV_overhead);
}

}