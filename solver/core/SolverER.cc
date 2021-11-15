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

inline void Solver::er_substitute(vec<Lit>& clause, struct ExtDefMap& extVarDefs) {
    // Get set of all literals in clause
    std::tr1::unordered_set<Lit> learntLits;
    for (int i = 1; i < clause.size(); i++) learntLits.insert(clause[i]);

    // Possible future investigation: should we ignore really long clauses in order to save time?
    // This would need a command-line option to toggle ignoring the clauses

    // Check for extension variables over literal pairs
    for (int i = 1; i < clause.size(); i++) {
        // Check if any extension variables are defined over this literal
        if (!extVarDefs.contains(clause[i])) continue;

        // TODO: would it be better to iterate over the extension definitions here?
        // Also consider sorting or doing a set intersection here
        for (int j = i + 1; j < clause.size(); j++) {
            // Check if any extension variables are defined over this literal pair
            std::tr1::unordered_map<std::pair<Lit, Lit>, Lit>::iterator it = extVarDefs.find(clause[i], clause[j]);
            if (it == extVarDefs.end()) continue;

            // Replace disjunction with intersection
            learntLits.erase(clause[i]);
            learntLits.erase(clause[j]);
            learntLits.insert(it->second);
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
    if (extVarDefs.size() > 0) {
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

// Add clause to extension definitions without performing simplification
void Solver::addExtDefClause(std::vector<CRef>& db, const vec<Lit>& ext_def_clause) {
    if (ext_def_clause.size() == 1){
        uncheckedEnqueue(ext_def_clause[0]);
    } else {
        CRef cr = ca.alloc(ext_def_clause, false);

        // Store clause in correct database
        db.push_back(cr);
        attachClause(cr);

        // Clauses containing extension variables should go in a separate database
        int numExtVarsInClause = getNumExtVars(ca[cr]);
        double extFrac = numExtVarsInClause / (double) ext_def_clause.size();
        extfrac_total += extFrac;

#if LBD_BASED_CLAUSE_DELETION
        const int clauseLBD = lbd(ca[cr]);
        ca[cr].activity() = clauseLBD;
        lbd_total += clauseLBD;
#endif
    }
}

inline std::vector<Var> Solver::er_add(
    std::tr1::unordered_map< Var, std::vector<CRef> >& er_def_db,
    struct ExtDefMap& er_def_map,
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
    vec<Lit> tmp;
    for (std::vector< std::pair< Var, std::pair<Lit, Lit> > >::const_iterator i = newDefMap.begin(); i != newDefMap.end(); i++) {
        // Get literals
        Lit x = mkLit(i->first);
        Lit a = i->second.first;
        Lit b = i->second.second;
        assert(var(x) > var(a) && var(x) > var(b));

        // Create extension clauses and add them to the extension definition database
        std::vector<CRef> defs;
        tmp.clear(); tmp.push(~x); tmp.push( a); tmp.push(b); addExtDefClause(defs, tmp);
        tmp.clear(); tmp.push( x); tmp.push(~a);              addExtDefClause(defs, tmp);
        tmp.clear(); tmp.push( x); tmp.push(~b);              addExtDefClause(defs, tmp);
        er_def_db.insert(std::make_pair(i->first, defs));

        // Propagate extension variable
        if (value(a) == l_False && value(b) == l_False) {
            uncheckedEnqueue(~x, defs[0]);
        } else if (value(a) == l_True) {
            uncheckedEnqueue(x, defs[1]);
        } else if (value(b) == l_True) {
            uncheckedEnqueue(x, defs[2]);
        }

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
    total_ext_vars += new_variables.size();
    extBuffer.clear();
    er_prioritize(new_variables);
    extTimerStop(ext_add_overhead);
}

static inline bool containsAny(Clause& c, const std::tr1::unordered_set<Var>& varSet) {
    for (int k = 0; k < c.size(); k++) {
        if (varSet.find(var(c[k])) != varSet.end()) {
            return true;
        }
    }
    return false;
}

std::tr1::unordered_set<Var> Solver::delExtVars(Minisat::vec<Minisat::CRef>& db, const std::tr1::unordered_set<Var>& varsToDeleteSet) {
    // Set of variables whose definitions cannot be deleted due to being locked
    std::tr1::unordered_set<Var> notDeleted;

    // Delete clauses which contain the extension variable
    // TODO: is there a more efficient way to implement this? e.g. have a list of clauses for each extension variable?
    // TODO: Should we queue ER clauses to be deleted when they become unlocked?
    int i, j;
    for (i = j = 0; i < db.size(); i++) {
        Clause& c = ca[db[i]];
        if (locked(c)) {
            db[j++] = db[i];
            // Add ext vars to notDeleted
            for (int k = 0; k < c.size(); k++) {
                if (varsToDeleteSet.find(var(c[k])) != varsToDeleteSet.end()) {
                    notDeleted.insert(var(c[k]));
                }
            }
        } else if (containsAny(c, varsToDeleteSet)) {
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
    return notDeleted;
}

std::tr1::unordered_set<Var> Solver::delExtVars(std::tr1::unordered_map< Var, std::vector<CRef> >& db, const std::tr1::unordered_set<Var>& varsToDelete) {
    std::tr1::unordered_set<Var> notDeleted;
    
    // Delete from variable definitions
    for (std::tr1::unordered_set<Var>::const_iterator i = varsToDelete.begin(); i != varsToDelete.end(); i++) {
        std::tr1::unordered_map< Var, std::vector<CRef> >::iterator it = extDefs.find(*i);
        std::vector<CRef>& defClauses = it->second;

        // Check if any of the clauses are locked
        bool canDelete = true;
        for (std::vector<CRef>::iterator k = defClauses.begin(); k != defClauses.end(); k++) {
            if (locked(ca[*k])) {
                canDelete = false;
                break;
            }
        }

        // Delete definition clauses
        if (canDelete) {
            for (std::vector<CRef>::iterator k = defClauses.begin(); k != defClauses.end(); k++) {
#if ER_USER_FILTER_HEURISTIC != ER_FILTER_HEURISTIC_NONE
                user_er_filter_delete_incremental(*k);
#endif
                removeClause(*k);
            }
        } else {
            notDeleted.insert(*i);
        }
    }
#if ER_USER_FILTER_HEURISTIC != ER_FILTER_HEURISTIC_NONE
    user_er_filter_delete_flush();
#endif

    return notDeleted;
}

static inline void setSubtract(std::tr1::unordered_set<Var>& a, const std::tr1::unordered_set<Var>& b) {
    for (std::tr1::unordered_set<Var>::const_iterator it = b.begin(); it != b.end(); it++) {
        a.erase(*it);
    }
}

void Solver::delExtVars(std::tr1::unordered_set<Var>(*er_delete_heuristic)(Solver&)) {
    // Delete variables
    extTimerStart();

    // Get variables to delete
    std::tr1::unordered_set<Var> varsToDelete = er_delete_heuristic(*this);

    // Option 1: delete all clauses containing the extension variables
    // Option 2: substitute extension variable with definition (TODO: unimplemented)

    // Only consider extension variables which are not part of the definition of another extension variable
    std::tr1::unordered_set<Var> cannotBeDeleted;
    for (std::tr1::unordered_set<Var>::iterator it = varsToDelete.begin(); it != varsToDelete.end(); it++)
        if (extVarDefs.contains(mkLit(*it, false)) || extVarDefs.contains(mkLit(*it, true)))
            cannotBeDeleted.insert(*it);
    setSubtract(varsToDelete, cannotBeDeleted);

    // Delete from learnt clauses
    cannotBeDeleted = delExtVars(extLearnts, varsToDelete);
    setSubtract(varsToDelete, cannotBeDeleted);

    // Delete clauses from extension definitions
    cannotBeDeleted = delExtVars(extDefs, varsToDelete);
    setSubtract(varsToDelete, cannotBeDeleted);

    // Remove variable definitions from other data structures
    deleted_ext_vars += varsToDelete.size();
    extVarDefs.erase(varsToDelete);

    extTimerStop(ext_delV_overhead);
}

}