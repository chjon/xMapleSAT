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

inline std::vector<Lit> Solver::er_substitute(vec<Lit>& clause, struct ExtDefMap& extVarDefs) {
    // Get indices of all literals that appear in an extension definition
    std::set<Lit> clauseLits;
    std::vector<int> defLitIndex;
    std::vector<Lit> substituted;
    for (int i = 0; i < clause.size(); i++) {
        clauseLits.insert(clause[i]);

        // Check if any extension variables are defined over this literal
        if (extVarDefs.contains(clause[i])) defLitIndex.push_back(i);
    }
    if (defLitIndex.size() <= 1) return substituted;

    // Check each pair of literals that appear in an extension definition
    int replaced = 0;
    for (unsigned int i = 0; i < defLitIndex.size(); i++) {
        if (clause[defLitIndex[i]] == lit_Undef) continue;
        for (unsigned int j = i + 1; j < defLitIndex.size(); j++) {
            if (clause[defLitIndex[j]] == lit_Undef) continue;

            // Check if any extension variables are defined over this literal pair
            std::tr1::unordered_map<std::pair<Lit, Lit>, Lit>::iterator it = extVarDefs.find(clause[defLitIndex[i]], clause[defLitIndex[j]]);
            if (it == extVarDefs.end()) continue;

            // Check if the extension variable already occurs in the clause
            std::set<Lit>::iterator it2 = clauseLits.find(it->second);
            if (it2 == clauseLits.end()) continue;

            // Replace the first literal with the extension literal and mark the second literal as invalid
            substituted.push_back(it->second);
            clause[defLitIndex[i]] = it->second;
            clause[defLitIndex[j]] = lit_Undef;
            replaced++;
            break;
        }
    }

    // Generate reduced learnt clause
    if (replaced > 0) {
        for (int i = 0, j = 0; i < clause.size(); i++) {
            if (clause[i] != lit_Undef) clause[j++] = clause[i];
        }
        clause.shrink(replaced);
    }
    return substituted;
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
            int clause_lbd = computeLBD(out_learnt);
            if (clause_lbd >= ext_min_lbd && clause_lbd <= ext_max_lbd)
#endif
            {
                // er_substitute returns a list of every extension literal that was substituted into the clause
                std::vector<Lit> extToPropagate = er_substitute(out_learnt, extVarDefs);

                // Extension variables might need to be propagated!
                // Their definition literals could be set without also setting the value of the extension variable
                // Propagate any extension variables that need to be propagated
                // This should be safe after backtracking
                for (unsigned int i = 0; i < extToPropagate.size(); i++) {
                    if (value(extToPropagate[i]) == l_Undef) {
                        // Check if any definition clauses are asserting
                        // We have to iterate through the clauses because we need the CRef for the reason clause when we propagate
                        std::vector<CRef> defClauses = extDefs.find(var(extToPropagate[i]))->second;
                        for (unsigned int j = 0; j < defClauses.size(); j++) {
                            Clause& c = ca[defClauses[j]];

                            // Find the extension variable and count the number of falsified non-extension literals
                            Lit extLit   = lit_Undef;
                            Lit falseLit = lit_Undef;
                            int numFalse = 0;
                            for (int k = 0; k < c.size(); k++) { // Definition clauses have size either 2 or 3
                                if (isExtVar(var(c[k]))) {
                                    extLit = c[k];
                                } else if (value(c[k]) == l_False) {
                                    falseLit = c[k];
                                    numFalse++;
                                }
                            }

                            // The clause will only be asserting if every literal but one is falsified
                            // This is an extension definition clause, so extLit should never be lit_Undef
                            assert(extLit != lit_Undef);
                            if (value(extLit) == l_Undef && numFalse == c.size() - 1) {
                                uncheckedEnqueue(extLit, level(var(falseLit)), defClauses[j]);
                                break;
                            }
                        }
                    }
                }

                // shiftUnassigned(out_learnt);

                // Shift unassigned literal to the beginning of the clause
                // There should be exactly one unassigned literal
                int i, j;
                for (i = j = 0; i < out_learnt.size(); i++) {
                    if (value(out_learnt[i]) == l_Undef) {
                        Lit tmp = out_learnt[i];
                        out_learnt[i] = out_learnt[j];
                        out_learnt[j] = tmp;
                        j++;
                    }
                }
                assert(j == 1);

                // Ensure the highest level is in index 1
                for (i = 1; i < out_learnt.size(); i++) {
                    if (level(var(out_learnt[i]) > level(var(out_learnt[1])))) {
                        Lit tmp = out_learnt[1];
                        out_learnt[1] = out_learnt[i];
                        out_learnt[i] = tmp;
                    }
                }
            }
        }
    }


    extTimerStop(ext_sub_overhead);
}

// Prioritize branching on a given set of literals
inline void Solver::er_prioritize(const std::vector<Var>& toPrioritize) {
    const double desiredActivityCHB   = activity_CHB  [order_heap_CHB[0]] * ext_prio_act;
    const double desiredActivityVSIDS = activity_VSIDS[order_heap_VSIDS[0]] * ext_prio_act;
    for (std::vector<Var>::const_iterator i = toPrioritize.begin(); i != toPrioritize.end(); i++) {
        Var v = *i;
        // Prioritize branching on our extension variables
        activity_CHB  [v] = desiredActivityCHB;
        activity_VSIDS[v] = desiredActivityVSIDS;

#if EXTENSION_FORCE_BRANCHING
        // This forces branching because of how branching is implemented when ANTI_EXPLORATION is turned on
        // FIXME: this only forces branching on the last extension variable we add here - maybe add a queue for force branch variables?
        canceled[v] = conflicts;
#endif
        if (order_heap_CHB  .inHeap(v)) order_heap_CHB  .decrease(v);
        if (order_heap_VSIDS.inHeap(v)) order_heap_VSIDS.decrease(v);
    }
}

void Solver::shiftUnassigned(vec<Lit>& clause) {
    if (clause.size() == 1) return;

    // Shift all unassigned literals to the beginning of the clause
    int i, j;
    for (i = j = 0; i < clause.size(); i++) {
        if (value(clause[i]) == l_Undef) {
            Lit tmp = clause[i]; clause[i] = clause[j]; clause[j] = tmp;
            j++;
        } else if (value(clause[i]) == l_True) {
            printf("Clause is satisfied!\n");
        }
    }

    const int num_unassigned = j;
    if (num_unassigned == 0) {
        // Ensure the highest level is in index 1 and the second-highest level is in index 0
        for (i = 0; i < clause.size(); i++) {
            if (level(var(clause[i])) >= level(var(clause[1]))) {
                Lit tmp = clause[0];
                clause[0] = clause[1];
                clause[1] = clause[i];
                clause[i] = tmp;
            }
        }
    } else if (num_unassigned == 1) {
        // Ensure the first literal has a value
        Lit tmp = clause[0];
        clause[0] = clause[1];
        clause[1] = tmp;
    }
}

// Add clause to extension definitions
// Assumes that the first literal is the new extension variable
int Solver::addExtDefClause(std::vector<CRef>& db, vec<Lit>& ext_def_clause) {
    // Try just using the default?
    // addClause_(ext_def_clause);

    if (ext_def_clause.size() == 1) {
        uncheckedEnqueue(ext_def_clause[0]);
    } else {
        // Make sure the first two literals are in the right order for the watchers
        shiftUnassigned(ext_def_clause);

        // Add clause to data structures
        CRef cr = ca.alloc(ext_def_clause, true);
        int lbd = computeLBD(ca[cr]);
        ca[cr].set_lbd(lbd);

        // Compute some stats
        int numExtVarsInClause = getNumExtVars(ca[cr]);
        double extFrac = numExtVarsInClause / (double) ext_def_clause.size();
        extfrac_total += extFrac;

        // Set initial clause LBD
        // lbd_total += lbd;

        // Add clause to db
        // clauses.push(cr);
        db.push_back(cr);

        attachClause(cr);

        // Do we ever actually need to propagate here?
        // BCP works by iterating through the literals on the trail 
        //
        // For ER_ADD_LOCATION_AFTER_RESTART:
        //    This means there are unit literals on the trail
        //    propagate() should handle this automatically
        //    
        // For ER_ADD_LOCATION_AFTER_CONFLICT:
        //    x = a v b: (-x a b)(x -a)(x -b)
        //    BCP will miss this if a and b were already set earlier
        //    We should backtrack to the appropriate level (max(lvl(a), lvl(b))) if we want to propagate, and
        //    let propagate() handle it for us

#if ER_USER_ADD_LOCATION == ER_ADD_LOCATION_AFTER_CONFLICT
        // Determine whether we need to propagate extension variable
        if (value(ext_def_clause[0]) == l_Undef) {
            int max_lvl = 0;
            bool allFalsified = true;
            for (int i = 1; i < ext_def_clause.size(); i++) {
                if (value(ext_def_clause[i]) != l_False) {
                    allFalsified = false;
                    break;
                }

                // Find the backtrack level in case we need to propagate
                int lvl = level(var(ext_def_clause[i]));
                if (lvl > max_lvl) max_lvl = lvl;
            }

            if (allFalsified) {
                return max_lvl;
            }
        }
#endif
    }
    return decisionLevel();
}

inline std::vector<Var> Solver::er_add(
    std::tr1::unordered_map< Var, std::vector<CRef> >& er_def_db,
    struct ExtDefMap& er_def_map,
    const std::vector< std::pair< Var, std::pair<Lit, Lit> > >& newDefMap,
    int& out_bt_level
) {
    // Add extension variables
    // It is the responsibility of the user heuristic to ensure that we do not have pre-existing extension variables
    // for the provided literal pairs

    // TODO: don't add the extension variable in the case where we have x1 = (a v b) and x2 = (x1 v a)
    // TODO: don't add the extension variable in the case where we have x1 = (a v b) and x2 = (x1 v -a)
    std::vector<Var> new_variables;
    for (unsigned int i = 0; i < newDefMap.size(); i++) {
        new_variables.push_back(newVar(ext_pref_sign));
    }

    // Add extension clauses
    int bt_level = decisionLevel();
    for (std::vector< std::pair< Var, std::pair<Lit, Lit> > >::const_iterator i = newDefMap.begin(); i != newDefMap.end(); i++) {
        // Get literals
        Lit x = mkLit(i->first);
        Lit a = i->second.first;
        Lit b = i->second.second;
        assert(var(x) > var(a) && var(x) > var(b));

        // Save definition
        er_def_map.insert(x, a, b);

        // Create extension clauses and add them to the extension definition database
        std::vector<CRef> defs;
        bt_level = std::min(bt_level, addExtDefClause(defs, ~x,  a,  b));
        bt_level = std::min(bt_level, addExtDefClause(defs,  x, ~a    ));
        bt_level = std::min(bt_level, addExtDefClause(defs,  x,     ~b));
        er_def_db.insert(std::make_pair(i->first, defs));
    }
    if (bt_level != decisionLevel()) out_bt_level = bt_level;
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
int Solver::addExtVars() {
    // Add the extension variables to our data structures
    int bt_level = -1;
    extTimerStart();
    const std::vector<Var> new_variables = er_add(extDefs, extVarDefs, extBuffer, bt_level);
    total_ext_vars += new_variables.size();
    extBuffer.clear();
    er_prioritize(new_variables);
    max_ext_vars = std::max(max_ext_vars, total_ext_vars - deleted_ext_vars);
    extTimerStop(ext_add_overhead);
    return bt_level;
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
    // cannotBeDeleted = delExtVars(extLearnts, varsToDelete);
    // setSubtract(varsToDelete, cannotBeDeleted);

    // Delete clauses from extension definitions
    cannotBeDeleted = delExtVars(extDefs, varsToDelete);
    setSubtract(varsToDelete, cannotBeDeleted);

    // Remove variable definitions from other data structures
    deleted_ext_vars += varsToDelete.size();
    extVarDefs.erase(varsToDelete);

    extTimerStop(ext_delV_overhead);
}

}