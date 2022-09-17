/***************************************************************************************[Solver.cc]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010, Niklas Sorensson
 
Chanseok Oh's MiniSat Patch Series -- Copyright (c) 2015, Chanseok Oh
 
Maple_LCM, Based on MapleCOMSPS_DRUP -- Copyright (c) 2017, Mao Luo, Chu-Min LI, Fan Xiao: implementing a learnt clause minimisation approach
Reference: M. Luo, C.-M. Li, F. Xiao, F. Manya, and Z. L. , “An effective learnt clause minimization approach for cdcl sat solvers,” in IJCAI-2017, 2017, pp. to–appear.

Maple_LCM_Dist, Based on Maple_LCM -- Copyright (c) 2017, Fan Xiao, Chu-Min LI, Mao Luo: using a new branching heuristic called Distance at the beginning of search

MapleLCMDistChronoBT, based on Maple_LCM_Dist -- Copyright (c), Alexander Nadel, Vadim Ryvchin: "Chronological Backtracking" in SAT-2018, pp. 111-121.

MapleLCMDistChronoBT-DL, based on MapleLCMDistChronoBT -- Copyright (c), Stepan Kochemazov, Oleg Zaikin, Victor Kondratiev, Alexander Semenov: The solver was augmented with heuristic that moves duplicate learnt clauses into the core/tier2 tiers depending on a number of parameters.

xMapleSAT -- Copyright (c) 2022, Jonathan Chung

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

#include <stdio.h>
#include "core/Solver.h"
#include "core/SolverER.h"
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

namespace Minisat {
    static const char* _ext = "EXT";
    static IntOption    opt_ext_freq       (_ext, "ext-freq","Number of conflicts to wait before trying to introduce an extension variable.\n", 2000, IntRange(0, INT32_MAX));
    static IntOption    opt_ext_wndw       (_ext, "ext-wndw","Number of clauses to consider when introducing extension variables.\n", 100, IntRange(0, INT32_MAX));
    static IntOption    opt_ext_num        (_ext, "ext-num", "Maximum number of extension variables to introduce at once\n", 1, IntRange(0, INT32_MAX));
    static DoubleOption opt_ext_prio       (_ext, "ext-prio","The fraction of maximum activity that should be given to new variables",  0.5, DoubleRange(0, false, HUGE_VAL, false));
    static BoolOption   opt_ext_sign       (_ext, "ext-sign","The default polarity of new extension variables (true = negative, false = positive)\n", true);
    static IntOption    opt_ext_min_width  (_ext, "ext-min-width", "Minimum clause width to select\n", 3, IntRange(0, INT32_MAX));
    static IntOption    opt_ext_max_width  (_ext, "ext-max-width", "Maximum clause width to select\n", 100, IntRange(0, INT32_MAX));
    static IntOption    opt_ext_sub_min_width  (_ext, "ext-sub-min-width", "Minimum width of clauses to substitute into\n", 3, IntRange(3, INT32_MAX));
    static IntOption    opt_ext_sub_max_width  (_ext, "ext-sub-max-width", "Maximum width of clauses to substitute into\n", INT32_MAX, IntRange(3, INT32_MAX));
    #if (ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_LBD) || ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LBD
    static IntOption    opt_ext_min_lbd    (_ext, "ext-min-lbd", "Minimum LBD of clauses to substitute into\n", 0, IntRange(0, INT32_MAX));
    static IntOption    opt_ext_max_lbd    (_ext, "ext-max-lbd", "Maximum LBD of clauses to substitute into\n", 5, IntRange(0, INT32_MAX));
    #endif
    #if ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY
    static DoubleOption opt_ext_act_thresh(_ext, "ext-act-thresh", "Activity threshold for extension variable deletion\n", 50, DoubleRange(1, false, HUGE_VAL, false));
    #elif ER_USER_DELETE_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY2
    static DoubleOption opt_ext_act_thresh(_ext, "ext-act-thresh", "Activity threshold for extension variable deletion\n", 0.5, DoubleRange(0, false, 1, false));
    #endif

    SolverER::SolverER(Solver* s)
        : total_ext_vars     (0)
        , deleted_ext_vars   (0)
        , max_ext_vars       (0)
        , conflict_extclauses(0)
        , learnt_extclauses  (0)
        , lbd_total          (0)
        , branchOnExt        (0)

        , prevExtensionConflict(0)

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
        , solver(s)
    {
        using namespace std::placeholders;

        // Bind filter heuristic
        #if ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_RANGE
            user_extFilPredicate = std::bind(&SolverER::user_extFilPredicate_width, this, _1);
        #elif ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LBD
            user_extFilPredicate = std::bind(&SolverER::user_extFilPredicate_lbd, this, _1);
        #endif

        // Bind clause selection heuristic
        #if ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_NONE
            user_extSelHeuristic = std::bind(&SolverER::user_extSelHeuristic_all, this, _1, _2, _3);
        #elif ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_ACTIVITY
            user_extSelHeuristic = std::bind(&SolverER::user_extSelHeuristic_activity, this, _1, _2, _3);
        #endif

        // Bind extension variable definition heuristic
        #if ER_USER_ADD_HEURISTIC == ER_ADD_HEURISTIC_SUBEXPR
            user_extDefHeuristic = std::bind(&SolverER::user_extDefHeuristic_subexpression, this, _1, _2, _3);
        #elif ER_USER_ADD_HEURISTIC == ER_ADD_HEURISTIC_RANDOM
            user_extDefHeuristic = std::bind(&SolverER::user_extDefHeuristic_random, this, _1, _2, _3);
        #endif

        // Bind clause substitution predicate
        user_extSubPredicate = std::bind(&SolverER::user_extSubPredicate_size_lbd, this, _1);

        // Bind variable deletion predicate
        #if ER_USER_DEL_HEURISTIC == ER_DELETE_HEURISTIC_NONE
            user_extDelPredicate = std::bind(&SolverER::user_extDelPredicate_none, this, _1);
        #elif ER_USER_DEL_HEURISTIC == ER_DELETE_HEURISTIC_ALL
            user_extDelPredicate = std::bind(&SolverER::user_extDelPredicate_all, this, _1);
        #elif ER_USER_DEL_HEURISTIC == ER_DELETE_HEURISTIC_ACTIVITY
            user_extDelPredicate = std::bind(&SolverER::user_extDelPredicate_activity, this, _1);
        #endif

        // Initialize overhead measurement
        initOverhead(ext_sel_overhead , 0, 0);
        initOverhead(ext_add_overhead , 0, 0);
        initOverhead(ext_delC_overhead, 0, 0);
        initOverhead(ext_delV_overhead, 0, 0);
        initOverhead(ext_sub_overhead , 0, 0);
        initOverhead(ext_stat_overhead, 0, 0);
    }

    SolverER::~SolverER() {}

    void SolverER::filterBatch(const vec<CRef>& candidates, FilterPredicate& filterPredicate) {
        extTimerStart();

        for (int i = 0; i < candidates.size(); i++) {
            CRef candidate = candidates[i];
            if (filterPredicate(candidate))
                m_filteredClauses.push_back(candidate);
        }

        extTimerStop(ext_sel_overhead);
    }

    void SolverER::filterIncremental(const CRef candidate, FilterPredicate& filterPredicate) {
        extTimerStart();
        
        if (filterPredicate(candidate))
            m_filteredClauses.push_back(candidate);
        
        extTimerStop(ext_sel_overhead);
    }

    void SolverER::selectClauses(SelectionHeuristic& selectionHeuristic) {
        filterBatch(solver->learnts_core , user_extFilPredicate);
        filterBatch(solver->learnts_tier2, user_extFilPredicate);

        extTimerStart();

        selectionHeuristic(m_selectedClauses, m_filteredClauses, ext_window);
        m_filteredClauses.clear();
        
        extTimerStop(ext_sel_overhead);
    }

    void SolverER::enforceWatcherInvariant(vec<Lit>& clause) const {
        assert(clause.size() > 1);

        // Swap all unassigned literals to the beginning of the clause
        int i, j;
        for (i = j = 0; i < clause.size(); i++) {
            if (value(clause[i]) == l_Undef) {
                Lit tmp = clause[i]; clause[i] = clause[j]; clause[j] = tmp;
                j++;
            }
        }

        const int num_unassigned = j;
        if (num_unassigned == 0) {
            // Move the two literals with the highest levels to the first two indices (O(n) partial selection sort)
            for (i = 0; i < 2; i++) {
                int maxLvl_j = i;
                for (j = i; j < clause.size(); j++) {
                    if (level(var(clause[j])) > level(var(clause[maxLvl_j]))) {
                        maxLvl_j = j;
                    }
                }
                Lit tmp = clause[i]; clause[i] = clause[maxLvl_j]; clause[maxLvl_j] = tmp;
            }

            // Swap the first two literals so that the highest level is in index 1 and the second-highest level is in index 0
            Lit tmp = clause[0]; clause[0] = clause[1]; clause[1] = tmp;

        } else if (num_unassigned == 1) {
            // Find the first literal assigned at the next-highest level:
            int maxLvl_j = 1;
            for (j = 2; j < clause.size(); j++)
                if (level(var(clause[j])) > level(var(clause[maxLvl_j])))
                    maxLvl_j = j;

            // Swap-in this literal at index 1:
            Lit tmp          = clause[maxLvl_j];
            clause[maxLvl_j] = clause[1];
            clause[1]        = tmp;
            // solver->cancelUntil(level(var(tmp)));
        }
    }

    void SolverER::defineExtVars(ExtDefHeuristic& extDefHeuristic) {
        extTimerStart();

        // Generate extension variable definitions
        extDefHeuristic(m_extVarDefBuffer, m_selectedClauses, ext_max_intro);

        // Update stats and clean up
        m_selectedClauses.clear();
        extTimerStop(ext_add_overhead);
    }

    void SolverER::introduceExtVars(std::tr1::unordered_map<Var, std::vector<CRef> >& ext_def_db) {
        if (m_extVarDefBuffer.size() == 0) return;

        extTimerStart();

        // Add extension variables
        // It is the responsibility of the user heuristic to ensure that we do not have pre-existing extension variables
        // for the provided literal pairs
        for (auto i = m_extVarDefBuffer.begin(); i != m_extVarDefBuffer.end(); i++) solver->newVar(ext_pref_sign);

        // Add extension definition clauses
        for (const ExtDef& def : m_extVarDefBuffer) {
            const Lit x = def.x, a = def.a, b = def.b;
            assert(var(x) >= originalNumVars && var(x) > var(a) && var(x) > var(b));

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
            ext_def_db.insert(std::make_pair(var(x), defs));
        }

        // TODO: Prioritize new variables
        // er_prioritize(m_extVarDefBuffer);

        // Update stats and clean up
        total_ext_vars += m_extVarDefBuffer.size();
        max_ext_vars = std::max(max_ext_vars, total_ext_vars - deleted_ext_vars);
        m_extVarDefBuffer.clear();

        extTimerStop(ext_add_overhead);
    }

    void SolverER::prioritize(const std::vector<ExtDef>& defs) {
        // FIXME: this only forces branching on the last extension variable we add here - maybe add a queue for force branch variables?
        const double desiredActivityCHB   = solver->activity_CHB  [solver->order_heap_CHB  [0]] * ext_prio_act;
        const double desiredActivityVSIDS = solver->activity_VSIDS[solver->order_heap_VSIDS[0]] * ext_prio_act;
        for (const ExtDef& def : defs) {
            Var v = var(def.x);
            // Prioritize branching on our extension variables
            solver->activity_CHB  [v] = desiredActivityCHB;
            solver->activity_VSIDS[v] = desiredActivityVSIDS;

#if EXTENSION_FORCE_BRANCHING
            // This forces branching because of how branching is implemented when ANTI_EXPLORATION is turned on
            solver->canceled[v] = conflicts;
#endif
            if (solver->order_heap_CHB  .inHeap(v)) solver->order_heap_CHB  .decrease(v);
            if (solver->order_heap_VSIDS.inHeap(v)) solver->order_heap_VSIDS.decrease(v);
        }
    }

    void SolverER::addExtDefClause(std::vector<CRef>& db, Lit ext_lit, vec<Lit>& ps) {
        assert(solver->decisionLevel() == 0);

        assert(solver->ok);

        // Check if clause is satisfied and remove false/duplicate literals:
        // TODO: make this optional depending on when we add the ext def clause
        sort(ps);
        Lit p; int i, j;

        if (solver->drup_file){
            solver->add_oc.clear();
            for (int i = 0; i < ps.size(); i++) solver->add_oc.push(ps[i]); }

        for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
            if (value(ps[i]) == l_True || ps[i] == ~p)
                return;
            else if (value(ps[i]) != l_False && ps[i] != p)
                ps[j++] = p = ps[i];
        ps.shrink(i - j);

        if (solver->drup_file && i != j){
    #ifdef BIN_DRUP
            solver->binDRUP('a', ps, solver->drup_file);
            solver->binDRUP('d', solver->add_oc, solver->drup_file);
    #else
            for (int i = 0; i < ps.size(); i++)
                fprintf(solver->drup_file, "%i ", (var(ps[i]) + 1) * (-2 * sign(ps[i]) + 1));
            fprintf(solver->drup_file, "0\n");

            fprintf(solver->drup_file, "d ");
            for (int i = 0; i < solver->add_oc.size(); i++)
                fprintf(solver->drup_file, "%i ", (var(solver->add_oc[i]) + 1) * (-2 * sign(solver->add_oc[i]) + 1));
            fprintf(solver->drup_file, "0\n");
    #endif
        }

        if (ps.size() == 0)
            return;
        else if (ps.size() == 1){
            solver->uncheckedEnqueue(ps[0]);
            solver->propagate();
            return;
        }else{
            CRef cr = solver->ca.alloc(ps, false);
            db.push_back(cr);
            solver->attachClause(cr);
        // }

        //////////////////////////////////////////////////////////////////////////////////

        // TODO: What happens if ER_USER_ADD_LOCATION == ER_ADD_LOCATION_AFTER_CONFLICT?
        // Do we need to propagate here?
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

        // if (ps.size() == 1) {
        //     solver->uncheckedEnqueue(ext_lit);
        // } else {
        //     // Make sure the first two literals are in the right order for the watchers
        //    enforceWatcherInvariant(ps);

        //     // Add clause to data structures
        //     ClauseAllocator& ca = solver->ca;
        //     CRef cr = ca.alloc(ps, false); // Allocating clause as if it were an original clause
        //     int lbd = solver->computeLBD(ca[cr]);
        //     ca[cr].set_lbd(lbd);

        //     // Add clause to db
        //     db.push_back(cr);
        //     solver->attachClause(cr);

            // Check whether the clause needs to be propagated
            if (value(ps[0]) == l_Undef && value(ps[1]) == l_False) {
                solver->uncheckedEnqueue(ps[0], level(var(ps[1])), cr);
            }
        }
    }

    void SolverER::getExtVarsToDelete(std::tr1::unordered_set<Lit>& varsToDelete, DeletionPredicate& deletionPredicate) const {
        // Iterate through current extension variables
        for (auto it = extDefs.begin(); it != extDefs.end(); it++) {
            Var x = it->first;

            // Find all extension variables that are not basis literals
            // Extension variables that participate in definitions of other variables cannot be deleted
            if (xdm.degree(mkLit(x, false)) > 0 || xdm.degree(mkLit(x, true)) > 0) continue;

            // Check whether the solver should delete the variable
            if (!deletionPredicate(x)) continue;

            // Check whether any of the extension definition clauses are not removable
            bool canDeleteDef = true;
            for (CRef cr : it->second) {
                Clause& c = solver->ca[cr];
                if (c.mark() == CORE || !c.removable() || solver->locked(c)) {
                    canDeleteDef = false;
                    break;
                }
            }

            // Queue extension variable to be deleted
            if (canDeleteDef) {
                varsToDelete.insert(mkLit(x));
            }
        }
    }

    void SolverER::deleteExtVars(DeletionPredicate& deletionPredicate) {
        extTimerStart();

        // Get extension variables to delete
        std::tr1::unordered_set<Lit> varsToDelete;
        getExtVarsToDelete(varsToDelete, deletionPredicate);

        // Delete clauses containing the extension variable
        // 1. TODO: Delete learnt clauses containing the extension variable
        //    NOTE: This is optional -- we can let the solver delete these by itself as clause activities decay 

        // 2. Delete extension variable definition clauses
        for (Lit x : varsToDelete) {
            for (CRef cr : extDefs.find(var(x))->second) {
                // TODO: also remove CRef from buffers (e.g. from incremental clause filtering)
                solver->removeClause(cr);
            }
            extDefs.erase(var(x));
        }

        // 3. Remove extension variable from definition map to prevent future clause substitution
        xdm.erase(varsToDelete); 

        // Update stats
        deleted_ext_vars += varsToDelete.size();
        extTimerStop(ext_delV_overhead);
    }

    void SolverER::relocAll(ClauseAllocator& to) {
        for (std::tr1::unordered_map< Var, std::vector<CRef> >::iterator it = extDefs.begin(); it != extDefs.end(); it++) {
            std::vector<CRef>& cs = it->second;
            unsigned int i, j;
            for (i = j = 0; i < cs.size(); i++) {
                // Reloc following example of clauses
                if (solver->ca[cs[i]].mark() != 1){
                    solver->ca.reloc(cs[i], to);
                    cs[j++] = cs[i]; }
            }
            cs.erase(cs.begin() + j, cs.end());
        }
    }
}