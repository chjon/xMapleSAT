/*******************************************************************************************[ERManager.cc]
xMaple*, extended resolution for Minisat-based solvers -- Copyright (c) 2022, Jonathan Chung, Vijay Ganesh, Sam Buss

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

#include "er/ERManager.h"

#include <stdio.h>
#include "core/Solver.h"
#include "mtl/Alg.h"
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

static inline void initOverhead(struct rusage& overhead, int s, int us) {
    overhead.ru_utime.tv_sec  = s;
    overhead.ru_utime.tv_usec = us;
}

using namespace Minisat;

static const char* _ext = "EXT";
static IntOption    opt_ext_freq       (_ext, "ext-freq","Number of conflicts to wait before trying to introduce an extension variable.\n", 2000, IntRange(0, INT32_MAX));
static IntOption    opt_ext_wndw       (_ext, "ext-wndw","Number of clauses to consider when introducing extension variables.\n", 100, IntRange(0, INT32_MAX));
static IntOption    opt_ext_num        (_ext, "ext-num", "Maximum number of extension variables to introduce at once\n", 1, IntRange(0, INT32_MAX));
static DoubleOption opt_ext_prio       (_ext, "ext-prio","The fraction of maximum activity that should be given to new variables",  1.0, DoubleRange(0, false, HUGE_VAL, false));
static BoolOption   opt_ext_sign       (_ext, "ext-sign","The default polarity of new extension variables (true = negative, false = positive)\n", true);
#if ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_RANGE
static IntOption    opt_ext_min_width  (_ext, "ext-min-width", "Minimum clause width to select\n", 3, IntRange(0, INT32_MAX));
static IntOption    opt_ext_max_width  (_ext, "ext-max-width", "Maximum clause width to select\n", 100, IntRange(0, INT32_MAX));
#endif
#if ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_WIDTH
static IntOption    opt_ext_sub_min_width  (_ext, "ext-sub-min-width", "Minimum width of clauses to substitute into\n", 3, IntRange(3, INT32_MAX));
static IntOption    opt_ext_sub_max_width  (_ext, "ext-sub-max-width", "Maximum width of clauses to substitute into\n", INT32_MAX, IntRange(3, INT32_MAX));
#endif
#if (ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_LBD) || ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LBD
static IntOption    opt_ext_min_lbd    (_ext, "ext-min-lbd", "Minimum LBD of clauses to substitute into\n", 0, IntRange(0, INT32_MAX));
static IntOption    opt_ext_max_lbd    (_ext, "ext-max-lbd", "Maximum LBD of clauses to substitute into\n", 5, IntRange(0, INT32_MAX));
#endif
#if ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY
static DoubleOption opt_ext_act_thresh(_ext, "ext-act-thresh", "Activity threshold for extension variable deletion\n", 50, DoubleRange(1, false, HUGE_VAL, false));
#elif ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY2
static DoubleOption opt_ext_act_thresh(_ext, "ext-act-thresh", "Activity threshold for extension variable deletion\n", 0.5, DoubleRange(0, false, 1, false));
#endif
static IntOption    opt_ext_del_freq(_ext, "ext-del-freq", "Number of conflicts to wait before trying to delete extension variables\n", 10000, IntRange(0, INT32_MAX));

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

ERManager::ERManager(Solver& s)
    ////////////////////
    // Solver references
    : assignmentTrail(s.assignmentTrail)
    , branchingHeuristicManager(s.branchingHeuristicManager)
    , clauseDatabase(s.clauseDatabase)
    , propagationQueue(s.propagationQueue)
    , randomNumberGenerator(s.randomNumberGenerator)
    , unitPropagator(s.unitPropagator)
    , ca(s.ca)
    , solver(s)

    ///////////////////
    // Member variables
    , prevExtensionConflict(0)
    , prevDelExtVarConflict(0)

    //////////////////////////
    // Command-line parameters
    , ext_freq         (opt_ext_freq)
    , ext_window       (opt_ext_wndw)
    , ext_max_intro    (opt_ext_num)
    , ext_prio_act     (opt_ext_prio)
    , ext_pref_sign    (opt_ext_sign)
#if ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_RANGE
    , ext_min_width    (opt_ext_min_width)
    , ext_max_width    (opt_ext_max_width)
#elif ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LONGEST
    , ext_filter_num   (opt_ext_filter_num)
#endif
#if ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_WIDTH
    , ext_sub_min_width (opt_ext_sub_min_width)
    , ext_sub_max_width (opt_ext_sub_max_width)
#endif
#if (ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_LBD) || ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LBD
    , ext_min_lbd      (opt_ext_min_lbd)
    , ext_max_lbd      (opt_ext_max_lbd)
#endif
#if ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY || ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY2
    , ext_act_threshold(opt_ext_act_thresh)
#endif
    , ext_del_freq(opt_ext_del_freq)

    /////////////
    // Statistics
    , total_ext_vars     (0)
    , tried_del_ext_vars (0)
    , deleted_ext_vars   (0)
    , max_ext_vars       (0)
    , conflict_extclauses(0)
    , learnt_extclauses  (0)
    , branchOnExt        (0)
{
    using namespace std::placeholders;

    // Bind filter heuristic
    #if ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_RANGE
        user_extFilPredicate = std::bind(&ERManager::user_extFilPredicate_width, this, _1);
    #elif ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LBD
        user_extFilPredicate = std::bind(&ERManager::user_extFilPredicate_lbd, this, _1);
    #endif

    // Bind clause selection heuristic
    #if ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_ALL
        user_extSelHeuristic = std::bind(&ERManager::user_extSelHeuristic_all, this, _1, _2, _3);
    #elif ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_ACTIVITY
        user_extSelHeuristic = std::bind(&ERManager::user_extSelHeuristic_activity, this, _1, _2, _3);
    #endif

    // Bind extension variable definition heuristic
    #if ER_USER_ADD_HEURISTIC == ER_ADD_HEURISTIC_SUBEXPR
        user_extDefHeuristic = std::bind(&ERManager::user_extDefHeuristic_subexpression, this, _1, _2, _3);
    #elif ER_USER_ADD_HEURISTIC == ER_ADD_HEURISTIC_RANDOM
        user_extDefHeuristic = std::bind(&ERManager::user_extDefHeuristic_random, this, _1, _2, _3);
    #endif

    // Bind clause substitution predicate
    user_extSubPredicate = std::bind(&ERManager::user_extSubPredicate_size_lbd, this, _1);

    // Bind variable deletion predicate setup function
    #if ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ALL
        user_extDelPredicateSetup = std::bind(&ERManager::user_extDelPredicateSetup_none, this);
    #elif ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY || ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY2
        user_extDelPredicateSetup = std::bind(&ERManager::user_extDelPredicateSetup_activity, this);
    #endif

    // Bind variable deletion predicate
    #if ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ALL
        user_extDelPredicate = std::bind(&ERManager::user_extDelPredicate_all, this, _1);
    #elif ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY || ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY2
        user_extDelPredicate = std::bind(&ERManager::user_extDelPredicate_activity, this, _1);
    #endif

    // Initialize overhead measurement
    initOverhead(ext_sel_overhead , 0, 0);
    initOverhead(ext_add_overhead , 0, 0);
    initOverhead(ext_delC_overhead, 0, 0);
    initOverhead(ext_delV_overhead, 0, 0);
    initOverhead(ext_sub_overhead , 0, 0);
    initOverhead(ext_stat_overhead, 0, 0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// CLAUSE SELECTION

void ERManager::filterIncremental(
    const CRef candidate,
    FilterPredicate& filterPredicate,
    HeuristicType heuristicType
) {
    extTimerStart();
    
    // Select data structure to use
    std::vector<CRef>* filteredClauses = nullptr;
    switch (heuristicType) {
        case HeuristicType::LER: filteredClauses = &m_filteredClauses_ler; break;
        default:                 filteredClauses = &m_filteredClauses    ; break;
    }

    // Add the clause if it satisfies the predicate
    if (filterPredicate(candidate))
        filteredClauses->push_back(candidate);
    
    extTimerStop(ext_sel_overhead);
}

void ERManager::selectClauses(
    SelectionHeuristic& selectionHeuristic,
    HeuristicType heuristicType, 
    unsigned int numKeepFiltered
) {
    extTimerStart();

    // Select data structure to use
    std::vector<CRef>* filteredClauses = nullptr;
    std::vector<CRef>* selectedClauses = nullptr;
    switch (heuristicType) {
        case HeuristicType::LER: {
            filteredClauses = &m_filteredClauses_ler;
            selectedClauses = &m_selectedClauses_ler;
        } break;
        default: {
            filteredClauses = &m_filteredClauses;
            selectedClauses = &m_selectedClauses;
        } break;
    }

    // For static metrics, it is preferable to use filterIncremental when learning clauses
    // instead of using selectLearntClauses here.
    // branchingHeuristicManager.selectLearntClauses(*filteredClauses, user_extFilPredicate);

    if (filteredClauses->size()) {
        selectionHeuristic(*selectedClauses, *filteredClauses, ext_window);
        const int end_index = std::max(
            0,
            static_cast<int>(filteredClauses->size()) - static_cast<int>(numKeepFiltered)
        );
        filteredClauses->erase(filteredClauses->begin(), filteredClauses->begin() + end_index);
    }
    extTimerStop(ext_sel_overhead);
}

////////////////////////////////////
// EXTENSION VARIABLE INTRODUCTION

void ERManager::defineExtVars(ExtDefHeuristic& extDefHeuristic, HeuristicType heuristicType) {
    extTimerStart();

    // Select data structure to use
    std::vector<CRef  >* selectedClauses = nullptr;
    std::vector<ExtDef>* extVarDefBuffer = nullptr;
    switch (heuristicType) {
        case HeuristicType::LER: {
            selectedClauses = &m_selectedClauses_ler;
            extVarDefBuffer = &m_extVarDefBuffer_ler;
        } break;
        default: {
            selectedClauses = &m_selectedClauses;
            extVarDefBuffer = &m_extVarDefBuffer;
        } break;
    }

    // Generate extension variable definitions
    extDefHeuristic(*extVarDefBuffer, *selectedClauses, ext_max_intro);

    // Update stats and clean up
    selectedClauses->clear();
    extTimerStop(ext_add_overhead);
}

void ERManager::introduceExtVars(
    std::tr1::unordered_map<Var, std::vector<CRef> >& ext_def_db,
    HeuristicType heuristicType
) {
    extTimerStart();

    // Select data structure to use
    std::vector<ExtDef>* extVarDefBuffer = nullptr;
    switch (heuristicType) {
        case HeuristicType::LER: {
            extVarDefBuffer = &m_extVarDefBuffer_ler;
        } break;
        default: {
            extVarDefBuffer = &m_extVarDefBuffer;
        } break;
    }

    // Add extension variables
    // It is the responsibility of the user heuristic to ensure that we do not have pre-existing extension variables
    // for the provided literal pairs
    // TODO: can we reuse the memory allocated for deleted variables? this should be safe as long as every
    // occurrence of the old variable has been removed
    for (auto i = extVarDefBuffer->begin(); i != extVarDefBuffer->end(); i++) solver.newVar(ext_pref_sign);

    // Add extension definition clauses
    for (const ExtDef& def : *extVarDefBuffer) {
        const Lit x = def.x, a = def.a, b = def.b;
        assert(var(x) >= originalNumVars && var(x) > var(a) && var(x) > var(b));
        
        // Update extension level
#if PRIORITIZE_ER || BUMP_ER
        extensionLevel[var(x)] = 1 + std::max(extensionLevel[var(a)], extensionLevel[var(b)]);
#endif

        // Save definition (x <=> a v b)
        xdm.insert(x, a, b);

        // Create extension clauses and save their IDs
        std::vector<CRef> defs;

        // Encode (x <=> a v b) as three clauses
        addExtDefClause(defs, x, {~x,  a,  b});
        addExtDefClause(defs, x, { x, ~a    });
        addExtDefClause(defs, x, { x,     ~b});

        // Introduce additional helper clauses
        for (const std::vector<Lit>& c : def.additionalClauses) addExtDefClause(defs, x, c);

        // Add extension clause IDs to the extension definition database
        extDefs.insert(std::make_pair(var(x), defs));
    }

    // Prioritize new variables
    // prioritize(*extVarDefBuffer);

    // Update stats and clean up
    total_ext_vars += extVarDefBuffer->size();
    max_ext_vars = std::max(max_ext_vars, total_ext_vars - deleted_ext_vars);
    extVarDefBuffer->clear();

    extTimerStop(ext_add_overhead);
}

void ERManager::prioritize(const std::vector<ExtDef>& defs) {
    // FIXME: this only forces branching on the last extension variable we add here - maybe add a
    // queue for force branch variables?

    // Find maximum activity
    double maxActivity = 0;
    const vec<double>& activities = branchingHeuristicManager.getActivityVSIDS();
    for (int i = 0; i < assignmentTrail.nVars(); i++)
        maxActivity = std::max(maxActivity, activities[i]);

    // Prioritize branching on our extension variables
    for (const ExtDef& def : defs) {
        Var v = var(def.x);
        branchingHeuristicManager.setActivity(v, maxActivity * ext_prio_act);
        branchingHeuristicManager.increasePriorityQueue(v);
    }
}

void ERManager::addExtDefClause(std::vector<CRef>& db, Lit ext_lit, vec<Lit>& ps) {
    // Copy clause
    sort(ps);
    Lit p; int i, j;
    for (i = j = 0, p = lit_Undef; i < ps.size(); i++) {
        // Don't add satisfied clauses
        if (
            ps[i] == ~p ||
            (assignmentTrail.value(ps[i]) == l_True && assignmentTrail.level(var(ps[i])) == 0)
        ) return;

        // Remove falsified literals
        if (
            ps[i] != p &&
            (assignmentTrail.value(ps[i]) != l_False || assignmentTrail.level(var(ps[i])) != 0)
        ) ps[j++] = p = ps[i];
    }
    ps.shrink(i - j);

    if (ps.size() == 0) {
        return;
    } else if (ps.size() == 1) {
        propagationQueue.enqueue(ps[0]);
        return;
    }

    // Enforce watcher invariant

    // Move undefined variables to the front
    for (i = j = 0; i < ps.size(); i++)
        if (assignmentTrail.value(ps[i]) == l_Undef)
            std::swap(ps[i], ps[j++]);

    // Move highest-level literal to ps[1] if there is only one unassigned variable
    for (i = j; j == 1 && i < ps.size(); i++)
        if (assignmentTrail.level(var(ps[i])) > assignmentTrail.level(var(ps[1])))
            std::swap(ps[i], ps[1]);

    // New clauses should not be conflicting!
    assert(assignmentTrail.value(ps[0]) != l_False);

    // Add clause
    CRef cr = clauseDatabase.addClause(ps, db, false);

    // Propagate clause if necessary
    if (assignmentTrail.value(ps[0]) == l_Undef && assignmentTrail.value(ps[1]) == l_False)
        propagationQueue.enqueue(ps[0], cr);
}

////////////////////////////////////
// EXTENSION VARIABLE SUBSTITUTION

void ERManager::substitute(vec<Lit>& clause, SubstitutionPredicate& p) {
    extTimerStart();
    vec<Lit>& extLits = tmp_vec; extLits.clear();

    // Nothing to do if we shouldn't substitute into the clause
    if (!p(clause)) goto EXIT_SUBSTITUTE;

    // Perform variable substitution
    xdm.substitute(clause, extLits);

    // Nothing to do if no variables were substituted
    if (extLits.size() == 0) goto EXIT_SUBSTITUTE;

    // Ensure variables are assigned so the clause is still asserting. Some extension variables
    // are undefined after substitution, so we need to propagate from their definitions
    for (int i = (extLits[0] == clause[0]) ? 1 : 0; i < extLits.size(); i++) {
        Lit x = extLits[i];
        if (assignmentTrail.value(x) != l_Undef) continue;

        int i_undef = -1, i_max = -1;
        CRef cr = findAssertingClause(i_undef, i_max, x, extDefs.find(var(x))->second);
        assert(cr != CRef_Undef);
        unitPropagator.enforceWatcherInvariant(cr, i_undef, i_max);
        Clause& c = ca[cr];

        propagationQueue.enqueue(c[0], cr);
    }

EXIT_SUBSTITUTE:;
    extTimerStop(ext_sub_overhead);
}

////////////////////////////////
// EXTENSION VARIABLE DELETION

void ERManager::getExtVarsToDelete(
    VarSet& varsToDelete,
    DeletionPredicate& deletionPredicate
) const {
    // Iterate through current extension variables
    for (auto it = extDefs.begin(); it != extDefs.end(); it++) {
        Var x = it->first;

        // Find all extension variables that are not basis literals
        // Extension variables that participate in definitions of other variables cannot be deleted
        if (xdm.degree(mkLit(x, false)) > 0 || xdm.degree(mkLit(x, true)) > 0) continue;

        // Check whether the solver should delete the variable
        if (!deletionPredicate(x)) continue;

        // Increment the number of attempted deletions
        tried_del_ext_vars++;

        // Check whether any of the extension definition clauses are not removable
        for (CRef cr : it->second) {
            Clause& c = ca[cr];
            if (assignmentTrail.locked(c))
                goto NEXT_VAR;
        }

        // Queue extension variable to be deleted
        varsToDelete.insert(x);

        NEXT_VAR:;
    }
}

void ERManager::deleteExtVars(DeletionPredicateSetup& setup, DeletionPredicate& deletionPredicate) {
    extTimerStart();

    // Get extension variables to delete
    VarSet varsToDelete;
    setup();
    getExtVarsToDelete(varsToDelete, deletionPredicate);

    extTimerStop(ext_delV_overhead);

    // Exit if there are no variables to delete
    if (varsToDelete.size() == 0) return;

    // Delete learnt clauses containing the extension variable
    // NOTE: This is optional -- we can let the solver delete these by itself as clause activities decay
    // extTimerStart();
    // deleteExtVarsFrom(solver->learnts, varsToDelete); // Delete from learnts
    // extTimerStop(ext_delC_overhead);

    extTimerStart();

    // Delete extension variable definition clauses
    for (Var x : varsToDelete) {
        for (CRef cr : extDefs.find(x)->second) {
            clauseDatabase.removeClause(cr);
        }
        extDefs.erase(x);
    }

    // Remove extension variable from definition map to prevent future clause substitution
    for (Var x : varsToDelete)
        xdm.erase(mkLit(x));

    // Flush clause selection buffers
    remove_flush();

    // Update stats
    deleted_ext_vars += varsToDelete.size();
    extTimerStop(ext_delV_overhead);
}

///////////////////////
// STATE MODIFICATION

void ERManager::relocAll(ClauseAllocator& to) {
    // Reloc CRefs stored in buffers
    for (unsigned int i = 0; i < m_filteredClauses.size(); i++) ca.reloc(m_filteredClauses[i], to);
    for (unsigned int i = 0; i < m_selectedClauses.size(); i++) ca.reloc(m_selectedClauses[i], to);
    for (unsigned int i = 0; i < m_filteredClauses_ler.size(); i++) ca.reloc(m_filteredClauses_ler[i], to);
    for (unsigned int i = 0; i < m_selectedClauses_ler.size(); i++) ca.reloc(m_selectedClauses_ler[i], to);

    // Reloc extension definition clauses
    for (std::tr1::unordered_map< Var, std::vector<CRef> >::iterator it = extDefs.begin(); it != extDefs.end(); it++) {
        std::vector<CRef>& cs = it->second;
        unsigned int i, j;
        for (i = j = 0; i < cs.size(); i++) {
            // Reloc following example of clauses
            if (ca[cs[i]].mark() != 1){
                ca.reloc(cs[i], to);
                cs[j++] = cs[i]; }
        }
        cs.erase(cs.begin() + j, cs.end());
    }
}

void ERManager::removeSatisfied() {
    // Iterate through every extension variable
    for (std::tr1::unordered_map< Var, std::vector<CRef> >::iterator it = extDefs.begin(); it != extDefs.end(); it++) {
        std::vector<CRef>& cs = it->second;

        // Iterate through the extension definition clauses
        unsigned int i, j;
        for (i = j = 0; i < cs.size(); i++) {
            Clause& c = ca[cs[i]];
            if (assignmentTrail.satisfied(c)) {
                clauseDatabase.removeClause(cs[i]);
            } else {
                cs[j++] = cs[i];
            }
        }

        // Shrink the vector if clauses were removed
        cs.erase(cs.begin() + j, cs.end());
    }

    remove_flush();
}

//////////////////////////////////////
// HELPER FUNCTIONS FOR SUBSTITUTION

CRef ERManager::findAssertingClause(int& i_undef, int& i_max, Lit x, std::vector<CRef>& cs) {
    // Find definition clause which asserts ~x
    int max_lvl;

    // Iterate through definition clauses
    for (CRef cr : cs) {
        Clause& c = ca[cr];
        i_undef = i_max = max_lvl = -1;

        // Check whether the clause is asserting
        for (int k = 0; k < c.size(); k++) {
            if (assignmentTrail.value(c[k]) == l_Undef) {
                if (c[k] != ~x) goto NEXT_CLAUSE;
                i_undef = k;
            } else if (assignmentTrail.value(c[k]) == l_True) {
                goto NEXT_CLAUSE;
            } else if (assignmentTrail.level(var(c[k])) > max_lvl) {
                max_lvl = assignmentTrail.level(var(c[k]));
                i_max = k;
            }
        }

        assert(i_undef != -1);
        return cr;
        NEXT_CLAUSE:;
    }

    // bug: the solver should never reach here if the extension definition clauses are present
    assert(false);
    return CRef_Undef;
}

/////////////////////////////////////////////////////
// HELPER FUNCTIONS FOR EXTENSION VARIABLE DELETION

/**
 * @brief Check whether a clause contains any variables in a set
 * 
 * @param c the clause to check
 * @param vars the set of variables
 * @return true iff the clause contains at least one variable which is in the set
 */
static inline bool containsAnyVar(Clause& c, VarSet& vars) {
    for (int i = 0; i < c.size(); i++)
        if (vars.find(var(c[i])) != vars.end())
            return true;

    return false;
}

void ERManager::deleteExtVarsFrom(vec<CRef>& db, VarSet& varsToDelete) {
    int i, j;
    for (i = j = 0; i < db.size(); i++) {
        CRef cr = db[i];
        Clause& c = ca[cr];
        if (!assignmentTrail.locked(c) && containsAnyVar(c, varsToDelete)) {
            clauseDatabase.removeClause(cr);
        } else {
            db[j++] = db[i];
        }
    }
    db.shrink(i - j);
}