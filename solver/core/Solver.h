/****************************************************************************************[Solver.h]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

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

#ifndef Minisat_Solver_h
#define Minisat_Solver_h

#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <vector>
#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "mtl/Alg.h"
#include "utils/System.h"
#include "utils/Options.h"
#include "core/SolverTypes.h"

#define MICROSEC_PER_SEC (1000000)

namespace Minisat {

struct LitPairMap {
    std::tr1::unordered_map< Lit, std::tr1::unordered_map<Lit, Lit> > map;
    std::tr1::unordered_map< Lit, std::pair<Lit, Lit> > xmap;

    inline bool contains (Lit a, Lit b) {
        std::tr1::unordered_map< Lit, std::tr1::unordered_map<Lit, Lit> >::const_iterator it1 = map.find(a);
        if (it1 == map.end()) return false;
        std::tr1::unordered_map<Lit, Lit>::const_iterator it2 = it1->second.find(b);
        return it2 != it1->second.end();
    }

    void insert (Lit x, Lit a, Lit b) {
        xmap.insert(std::make_pair(x, std::make_pair(a, b)));

        // Insert for tuple <a, b>
        std::tr1::unordered_map< Lit, std::tr1::unordered_map<Lit, Lit> >::iterator it1 = map.find(a);
        if (it1 == map.end()) {
            std::tr1::unordered_map<Lit, Lit> submap;
            submap.insert(std::make_pair(b, x));
            map.insert(std::make_pair(a, submap));
        } else {
            it1->second.insert(std::make_pair(b, x));
        }

        // Insert for tuple <b, a>
        std::tr1::unordered_map< Lit, std::tr1::unordered_map<Lit, Lit> >::iterator it2 = map.find(b);
        if (it2 == map.end()) {
            std::tr1::unordered_map<Lit, Lit> submap;
            submap.insert(std::make_pair(a, x));
            map.insert(std::make_pair(b, submap));
        } else {
            it2->second.insert(std::make_pair(a, x));
        }
    }

    inline void erase(Lit a, Lit b) {
        // Check if the pair is in the map
        std::tr1::unordered_map< Lit, std::tr1::unordered_map<Lit, Lit> >::iterator it1 = map.find(a);
        if (it1 == map.end()) return;
        std::tr1::unordered_map<Lit, Lit>& submap = it1->second;
        std::tr1::unordered_map<Lit, Lit>::const_iterator it2 = submap.find(b);
        if (it2 == submap.end()) return;

        // Delete the pair
        if (submap.size() == 0) map.erase(it1);
        else                    submap.erase(it2);

        // Delete the pair in the reverse order
        it1 = map.find(b);
        submap = it1->second;
        it2 = submap.find(a);
        if (submap.size() == 0) map.erase(it1);
        else                    submap.erase(it2);
    }

    void erase(const std::tr1::unordered_set<Var>& defsToDelete) {
        for (std::tr1::unordered_set<Var>::const_iterator i = defsToDelete.begin(); i != defsToDelete.end(); i++) {
            std::tr1::unordered_map< Lit, std::pair<Lit, Lit> >::iterator it = xmap.find(mkLit(*i));
            std::pair<Lit, Lit>& def = it->second;
            erase(def.first, def.second);
        }
    }
};

//=================================================================================================
// Solver -- the main class:

class Solver {
public:

    // Constructor/Destructor:
    //
    Solver();
    virtual ~Solver();

    // Problem specification:
    //
    Var     newVar    (bool polarity = true, bool dvar = true); // Add a new variable with parameters specifying variable mode.

    bool    addClause (const vec<Lit>& ps);                     // Add a clause to the solver. 
    bool    addEmptyClause();                                   // Add the empty clause, making the solver contradictory.
    bool    addClause (Lit p);                                  // Add a unit clause to the solver. 
    bool    addClause (Lit p, Lit q);                           // Add a binary clause to the solver. 
    bool    addClause (Lit p, Lit q, Lit r);                    // Add a ternary clause to the solver. 
    bool    addClause_(      vec<Lit>& ps);                     // Add a clause to the solver without making superflous internal copy. Will
                                                                // change the passed vector 'ps'.
    bool    addClauseToDB(vec<CRef>& db, Lit p);                // Add a unit clause to the solver. 
    bool    addClauseToDB(vec<CRef>& db, Lit p, Lit q);         // Add a binary clause to the solver. 
    bool    addClauseToDB(vec<CRef>& db, Lit p, Lit q, Lit r);  // Add a ternary clause to the solver. 
    bool    addClauseToDB(vec<CRef>& clauseDB, vec<Lit>& ps);   // Add a clause to a specific clause DB without making superflous internal
                                                                // copy. Will change the passed vector 'ps'.

    // Solving:
    //
    bool    simplify     ();                        // Removes already satisfied clauses.
    bool    solve        (const vec<Lit>& assumps); // Search for a model that respects a given set of assumptions.
    lbool   solveLimited (const vec<Lit>& assumps); // Search for a model that respects a given set of assumptions (With resource constraints).
    bool    solve        ();                        // Search without assumptions.
    bool    solve        (Lit p);                   // Search for a model that respects a single assumption.
    bool    solve        (Lit p, Lit q);            // Search for a model that respects two assumptions.
    bool    solve        (Lit p, Lit q, Lit r);     // Search for a model that respects three assumptions.
    bool    okay         () const;                  // FALSE means solver is in a conflicting state

    void    toDimacs     (FILE* f, const vec<Lit>& assumps);            // Write CNF to file in DIMACS-format.
    void    toDimacs     (const char *file, const vec<Lit>& assumps);
    void    toDimacs     (FILE* f, Clause& c, vec<Var>& map, Var& max);

    // Convenience versions of 'toDimacs()':
    void    toDimacs     (const char* file);
    void    toDimacs     (const char* file, Lit p);
    void    toDimacs     (const char* file, Lit p, Lit q);
    void    toDimacs     (const char* file, Lit p, Lit q, Lit r);
    
    // Variable mode:
    // 
    void    setPolarity    (Var v, bool b); // Declare which polarity the decision heuristic should use for a variable. Requires mode 'polarity_user'.
    void    setDecisionVar (Var v, bool b); // Declare if a variable should be eligible for selection in the decision heuristic.

    // Read state:
    //
    lbool   value      (Var x) const;       // The current value of a variable.
    lbool   value      (Lit p) const;       // The current value of a literal.
    lbool   modelValue (Var x) const;       // The value of a variable in the last model. The last call to solve must have been satisfiable.
    lbool   modelValue (Lit p) const;       // The value of a literal in the last model. The last call to solve must have been satisfiable.
    int     nAssigns   ()      const;       // The current number of assigned literals.
    int     nClauses   ()      const;       // The current number of original clauses.
    int     nLearnts   ()      const;       // The current number of learnt clauses.
    int     nVars      ()      const;       // The current number of variables.
    int     nFreeVars  ()      const;

    // EXTENDED RESOLUTION - statistics
    int     nExtLearnts()      const;       // The current number of learnt extension clauses.
    int     nExtDefs   ()      const;       // The current number of extension definition clauses.
    int     nExtVars   ()      const;       // The current number of extension variables.

    // Resource contraints:
    //
    void    setConfBudget(int64_t x);
    void    setPropBudget(int64_t x);
    void    budgetOff();
    void    interrupt();          // Trigger a (potentially asynchronous) interruption of the solver.
    void    clearInterrupt();     // Clear interrupt indicator flag.
    bool    interrupted() const; // Check if the solver has been interrupted

    // Memory managment:
    //
    virtual void garbageCollect();
    void    checkGarbage(double gf);
    void    checkGarbage();

    // Extra results: (read-only member variable)
    //
    vec<lbool> model;             // If problem is satisfiable, this vector contains the model (if any).
    vec<Lit>   conflict;          // If problem is unsatisfiable (possibly under assumptions),
                                  // this vector represent the final conflict clause expressed in the assumptions.

    // Mode of operation:
    //
    int       verbosity;
#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
    double    step_size;
    double    step_size_dec;
    double    min_step_size;
#endif
#if BRANCHING_HEURISTIC == VSIDS
    double    var_decay;
#endif
#if ! LBD_BASED_CLAUSE_DELETION
    double    clause_decay;
#endif
    double    random_var_freq;
    double    random_seed;
    bool      luby_restart;
    int       ccmin_mode;         // Controls conflict clause minimization (0=none, 1=basic, 2=deep).
    int       phase_saving;       // Controls the level of phase saving (0=none, 1=limited, 2=full).
    bool      rnd_pol;            // Use random polarities for branching heuristics.
    bool      rnd_init_act;       // Initialize variable activities with a small random value.
    double    garbage_frac;       // The fraction of wasted memory allowed before a garbage collection is triggered.

    int       restart_first;      // The initial restart limit.                                                                (default 100)
    double    restart_inc;        // The factor with which the restart limit is multiplied in each restart.                    (default 1.5)

    int       ext_freq;           // Number of conflicts to wait before trying to introduce an extension variable              (default 2000)
    int       ext_window;         // Number of clauses to consider when introducing extension variables.                       (default 100)
    int       ext_max_intro;      // Maximum number of extension variables to introduce at once.                               (default 1)
#if ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_RANGE
    int       ext_min_width;      // Minimum clause width to consider when selecting clauses
    int       ext_max_width;      // Maximum clause width to consider when selecting clauses
#elif ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LONGEST
    int       ext_filter_num;     // Maximum number of clauses after the filter step
#endif
#if ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_WIDTH
    int       ext_sub_min_width;      // Minimum width of clauses to substitute into
    int       ext_sub_max_width;      // Maximum width of clauses to substitute into
#endif
#if (ER_USER_SUBSTITUTE_HEURISTIC & ER_SUBSTITUTE_HEURISTIC_LBD) || ER_USER_FILTER_HEURISTIC == ER_FILTER_HEURISTIC_LBD
    int       ext_min_lbd;        // Minimum LBD of clauses to substitute into
    int       ext_max_lbd;        // Maximum LBD of clauses to substitute into
#endif
    double    learntsize_factor;  // The intitial limit for learnt clauses is a factor of the original clauses.                (default 1 / 3)
    double    learntsize_inc;     // The limit for learnt clauses is multiplied with this factor each restart.                 (default 1.1)

    int       learntsize_adjust_start_confl;
    double    learntsize_adjust_inc;

    // Statistics: (read-only member variable)
    //
    uint64_t solves, starts, decisions, rnd_decisions, propagations, conflicts;
    uint64_t dec_vars, clauses_literals, learnts_literals, max_literals, tot_literals;
    
    // EXTENDED RESOLUTION - statistics
    // read-only member variables
    uint64_t conflict_extclauses, learnt_extclauses, lbd_total, branchOnExt;
    double extfrac_total;
    struct rusage ext_timer_start, ext_timer_end;
    struct rusage ext_sel_overhead; // Overhead for selecting clauses for adding extension variables
    struct rusage ext_add_overhead; // Overhead for adding extension variables
    struct rusage ext_delC_overhead; // Overhead for deleting clauses containing extension variables
    struct rusage ext_delV_overhead; // Overhead for deleting extension variables
    struct rusage ext_sub_overhead; // Overhead for substituting disjunctions containing extension variables
    struct rusage ext_stat_overhead; // Overhead for measuring statistics
    double extTimerRead(unsigned int i); // 0: sel, 1: add, 2: delC, 3: delV, 4: sub, 5: stat

    uint64_t lbd_calls;
    vec<uint64_t> lbd_seen;
    vec<uint64_t> picked;
    vec<uint64_t> conflicted;
#if ALMOST_CONFLICT
    vec<uint64_t> almost_conflicted;
#endif
#if ANTI_EXPLORATION
    vec<uint64_t> canceled;
#endif
#if BRANCHING_HEURISTIC == CHB
    vec<uint64_t> last_conflict;
    int action;
    double reward_multiplier;
#endif

    vec<long double> total_actual_rewards;
    vec<int> total_actual_count;

protected:

    // Helper structures:
    //
    struct VarData { CRef reason; int level; };
    static inline VarData mkVarData(CRef cr, int l){ VarData d = {cr, l}; return d; }

    struct Watcher {
        CRef cref;
        Lit  blocker;
        Watcher(CRef cr, Lit p) : cref(cr), blocker(p) {}
        bool operator==(const Watcher& w) const { return cref == w.cref; }
        bool operator!=(const Watcher& w) const { return cref != w.cref; }
    };

    struct WatcherDeleted
    {
        const ClauseAllocator& ca;
        WatcherDeleted(const ClauseAllocator& _ca) : ca(_ca) {}
        bool operator()(const Watcher& w) const { return ca[w.cref].mark() == 1; }
    };

    struct VarOrderLt {
        const vec<double>&  activity;
        bool operator () (Var x, Var y) const { return activity[x] > activity[y]; }
        VarOrderLt(const vec<double>&  act) : activity(act) { }
    };

    // Solver state:
    //
    bool                ok;               // If FALSE, the constraints are already unsatisfiable. No part of the solver state may be used!
    vec<CRef>           clauses;          // List of problem clauses.
    vec<CRef>           learnts;          // List of learnt clauses.

    // EXTENDED RESOLUTION - clause databases
    vec<CRef>           extLearnts;       // List of learnt extension clauses (learnt clauses which contain extension variables).
    vec<CRef>           extDefs;          // List of extension definition clauses.

    std::vector< std::pair< Var, std::pair<Lit, Lit> > > extBuffer; // Buffer of extension variable definitions to add

#if ! LBD_BASED_CLAUSE_DELETION
    double              cla_inc;          // Amount to bump next clause with.
#endif
    vec<double>         activity;         // A heuristic measurement of the activity of a variable.
    double              var_inc;          // Amount to bump next variable with.
    OccLists<Lit, vec<Watcher>, WatcherDeleted>
                        watches;          // 'watches[lit]' is a list of constraints watching 'lit' (will go there if literal becomes true).
    vec<lbool>          assigns;          // The current assignments.
    vec<char>           polarity;         // The preferred polarity of each variable.
    vec<char>           decision;         // Declares if a variable is eligible for selection in the decision heuristic.
    vec<Lit>            trail;            // Assignment stack; stores all assigments made in the order they were made.
    vec<int>            trail_lim;        // Separator indices for different decision levels in 'trail'.
    vec<VarData>        vardata;          // Stores reason and level for each variable.
    int                 qhead;            // Head of queue (as index into the trail -- no more explicit propagation queue in MiniSat).
    int                 simpDB_assigns;   // Number of top-level assignments since last execution of 'simplify()'.
    int64_t             simpDB_props;     // Remaining number of propagations that must be made before next execution of 'simplify()'.
    vec<Lit>            assumptions;      // Current set of assumptions provided to solve by the user.
    Heap<VarOrderLt>    order_heap;       // A priority queue of variables ordered with respect to the variable activity.
    double              progress_estimate;// Set by 'search()'.
    bool                remove_satisfied; // Indicates whether possibly inefficient linear scan for satisfied clauses should be performed in 'simplify'.
    
    // EXTENDED RESOLUTION - solver state
    struct LitPairMap extVarDefs;
                                             // Extension variable definitions - key is a pair of literals and value is the corresponding extension variable
                                             // This map is used for replacing disjunctions with the corresponding extension variable
                                             // This is NOT the same as the extension variable introduction heuristic
#if ER_USER_FILTER_HEURISTIC != ER_FILTER_HEURISTIC_NONE
    std::tr1::unordered_set<CRef> er_deletedClauses;
    std::vector<CRef> er_filteredClauses;
                                             // List of clauses which can be selected by the clause selection heuristic
                                             // This represents the result of an initial filtering step, such as filtering by clause width
                                             // Special care needs to be taken while deleting clauses
#endif
    int               originalNumVars;       // The number of variables in the original formula
                                             // This value is used to quickly check whether a variable is an extension variable
    long unsigned int prevExtensionConflict; // Stores the last time extension variables were added
                                             // This is used to check whether to run the extension variable introduction heuristic after a restart

    ClauseAllocator     ca;

    // Temporaries (to reduce allocation overhead). Each variable is prefixed by the method in which it is
    // used, except 'seen' which is used in several places.
    //
    vec<char>           seen;
    vec<Lit>            analyze_stack;
    vec<Lit>            analyze_toclear;
    vec<Lit>            add_tmp;
    std::vector<Lit>    make_tmp;

    double              max_learnts;
    double              learntsize_adjust_confl;
    int                 learntsize_adjust_cnt;

    // Resource contraints:
    //
    int64_t             conflict_budget;    // -1 means no budget.
    int64_t             propagation_budget; // -1 means no budget.
    bool                asynch_interrupt;

    // Main internal methods:
    //
    void     insertVarOrder   (Var x);                                                 // Insert a variable in the decision order priority queue.
    Lit      pickBranchLit    ();                                                      // Return the next decision variable.
    void     newDecisionLevel ();                                                      // Begins a new decision level.
    void     uncheckedEnqueue (Lit p, CRef from = CRef_Undef);                         // Enqueue a literal. Assumes value of literal is undefined.
    bool     enqueue          (Lit p, CRef from = CRef_Undef);                         // Test if fact 'p' contradicts current state, enqueue otherwise.
    CRef     propagate        ();                                                      // Perform unit propagation. Returns possibly conflicting clause.
    void     cancelUntil      (int level);                                             // Backtrack until a certain level.
    void     analyze          (CRef confl, vec<Lit>& out_learnt, int& out_btlevel);    // (bt = backtrack)
    void     analyzeFinal     (Lit p, vec<Lit>& out_conflict);                         // COULD THIS BE IMPLEMENTED BY THE ORDINARIY "analyze" BY SOME REASONABLE GENERALIZATION?
    bool     litRedundant     (Lit p, uint32_t abstract_levels);                       // (helper method for 'analyze()')

    template<class V> int lbd (const V& clause) {
        lbd_calls++;
        int lbd = 0;
        for (int i = 0; i < clause.size(); i++) {
            int l = level(var(clause[i]));
            if (lbd_seen[l] != lbd_calls) {
                lbd++;
                lbd_seen[l] = lbd_calls;
            }
        }
        return lbd;
    }
    lbool    search           (int nof_conflicts);                                     // Search for a given number of conflicts.
    lbool    solve_           ();                                                      // Main solve method (assumptions given in 'assumptions').
    void     reduceDB         ();                                                      // Reduce the set of learnt clauses.
    void     reduceDB         (Minisat::vec<Minisat::CRef>& db);                       // Reduce the provided set of learnt clauses.
    void     removeSatisfied  (vec<CRef>& cs);                                         // Shrink 'cs' to contain only non-satisfied clauses.
    void     rebuildOrderHeap ();

    // Maintaining Variable/Clause activity:
    //
#if BRANCHING_HEURISTIC == VSIDS
    void     varDecayActivity ();                      // Decay all variables with the specified factor. Implemented by increasing the 'bump' value instead.
    void     varBumpActivity  (Var v, double inc);     // Increase a variable with the current 'bump' value.
    void     varBumpActivity  (Var v);                 // Increase a variable with the current 'bump' value.
#endif
#if ! LBD_BASED_CLAUSE_DELETION
    void     claDecayActivity ();                      // Decay all clauses with the specified factor. Implemented by increasing the 'bump' value instead.
    void     claBumpActivity  (Clause& c);             // Increase a clause with the current 'bump' value.
#endif

    // Operations on clauses:
    //
    void     attachClause     (CRef cr);               // Attach a clause to watcher lists.
    void     detachClause     (CRef cr, bool strict = false); // Detach a clause to watcher lists.
    void     removeClause     (CRef cr);               // Detach and free a clause.
    bool     locked           (const Clause& c) const; // Returns TRUE if a clause is a reason for some implication in the current state.
    bool     satisfied        (const Clause& c) const; // Returns TRUE if a clause is satisfied in the current state.

    void     relocHelper      (CRef& cr, ClauseAllocator& to, std::tr1::unordered_set<CRef>& newFilteredClauses);
    void     relocAll         (ClauseAllocator& to);

    // Misc:
    //
    int      decisionLevel    ()      const; // Gives the current decisionlevel.
    uint32_t abstractLevel    (Var x) const; // Used to represent an abstraction of sets of decision levels.
    CRef     reason           (Var x) const;
    int      level            (Var x) const;
    double   progressEstimate ()      const; // DELETE THIS ?? IT'S NOT VERY USEFUL ...
    bool     withinBudget     ()      const;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // EXTENDED RESOLUTION - internal solver functions
    // check the type of clause or variable
    int getNumExtVars (const Clause& c) const; // Count the number of extension variables in a clause
    int isExtVar      (Var x)           const; // Check whether a variable is an extension variable - returns 1 if it is and 0 if it is not
    
    // Main internal methods

    // Description:
    //   Pick extension variable definitions and add them to the extension definition buffer.
    //
    // Parameters:
    //   er_select_heuristic:
    //     A heuristic for picking candidate clauses to be used for the extended resolution variable
    //     introduction heuristic.
    //   er_add_heuristic:
    //     A heuristic for identifying extension variable definitions.
    //   numClausesToConsider:
    //     The number of clauses to consider when looking for new variable definitions.
    //   maxNumNewVars:
    //     the maximum number of new extension variables to introduce.
    void generateExtVars (
        std::vector<CRef>(*er_select_heuristic)(Solver&, unsigned int),
        std::vector< std::pair< Var, std::pair<Lit, Lit> > >(*er_add_heuristic)(Solver&, std::vector<CRef>&, unsigned int),
        unsigned int numClausesToConsider,
        unsigned int maxNumNewVars
    );

    // Description:
    //   Add extension variables from the extension definition buffer to our data structures and prioritize branching on them.
    void addExtVars ();

    // Internal helpers for addExtVars
    void er_prioritize(const std::vector<Var>& toPrioritize);
    std::vector<Var> er_add(
        vec<CRef>& er_def_db,
        struct LitPairMap& er_def_map,
        const std::vector< std::pair< Var, std::pair<Lit, Lit> > >& newDefMap
    );

    // Description:
    //   Select extension variables to delete and delete them
    //
    // Parameters:
    //   er_delete_heuristic:
    //     A heuristic for identifying extension variables to delete.
    void delExtVars (std::tr1::unordered_set<Var>(*er_delete_heuristic)(Solver&));

    // Description:
    //   Delete a specified set of extension variables from a given clause database. This is a helper function
    //   for the overloaded delExtVars function above
    //
    // Parameters:
    //   db     : The clause database to delete from
    //   extvars: The set of extension variables to delete
    void delExtVars (Minisat::vec<Minisat::CRef>& db, const std::tr1::unordered_set<Var>& extvars);

    // Description:
    //   Replace variable disjunctions in candidate learnt clauses with the corresponding extension variable
    //   if one exists. e.g. if x = a v b, and we learn a clause C = (a, b, c, d, ...), replace a, b with x and
    //   learn the clause C' = (x, c, d, ...) instead.
    //
    // Parameters:
    //   out_learnt: The candidate clause to be learnt
    void substituteExt (vec<Lit>& out_learnt);

    // Internal helper for substituteExt
    static void er_substitute(
        vec<Lit>& out_learnt,
        struct LitPairMap& extVarDefs
    );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // EXTENDED RESOLUTION - user functions/heuristics
    // For organization, all these function names should be prefixed according to their usage.
    //
    // Prefixes:
    //   user_er_select_
    //   user_er_add_
    //   user_er_delete_

    // FIXME: We want the solver object to be const when passed into the heuristic functions. Passing in a mutable solver and
    //   inspecting its internal state from the heuristic functions is bad API design. Currently, the solver needs to be
    //   mutable because we need access to its clause allocator when dereferencing CRef pointers, and the clause allocator
    //   index operator returns mutable state. Can we overload the operator with a const version? Investigate this.

    ///// [ user_er_select_ ] /////
    // Description:
    //   User functions for picking candidate clauses for the extended resolution variable introduction heuristic.
    //
    // Parameters:
    //   solver    : the solver to select clauses from
    //   numClauses: the number of clauses to select
    //
    // Return value:
    //   The function should return a list of clauses which are somehow interesting for ER.

    // Activity-based clause selection - select the most active clauses. This heuristic focuses on locality.
    static void copy_k_largest_activity(vec<CRef>& target, vec<CRef>& source, Solver& solver, unsigned int k);

    // Quickselect based on clause activity
    static int  partition_activity(vec<CRef>& db, ClauseAllocator& ca, int l, int r, int pivot);
    static void quickselect_activity(vec<CRef>& db, Solver& solver, int l, int r, int k);
    static int  partition_count(std::vector< std::pair<Lit, Lit> >& db, std::tr1::unordered_map<std::pair<Lit, Lit>, int>& subexpr_count, int l, int r, int pivot);
    static void quickselect_count(std::vector< std::pair<Lit, Lit> >& db, std::tr1::unordered_map<std::pair<Lit, Lit>, int>& subexpr_count, Solver& solver, int l, int r, int k);

    void user_er_filter_incremental(const CRef candidate);
    static void user_er_select_filter_widths(vec<CRef>& output, const vec<CRef>& clauses, ClauseAllocator& ca, int minWidth, int maxWidth);
    void user_er_filter_delete_incremental(CRef cr);
    void user_er_filter_delete_flush(void);

    void user_er_reloc(ClauseAllocator& to);

    static std::vector<CRef> user_er_select          (Solver& solver, unsigned int numClauses);
#if ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_NONE
    static std::vector<CRef> user_er_select_naive    (Solver& solver, unsigned int numClauses);
#elif ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_ACTIVITY
    static std::vector<CRef> user_er_select_activity (Solver& solver, unsigned int numClauses);
#elif ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_ACTIVITY2
    static std::vector<CRef> user_er_select_activity2(Solver& solver, unsigned int numClauses);
#elif ER_USER_SELECT_HEURISTIC == ER_SELECT_HEURISTIC_GLUCOSER
    static std::vector<CRef> user_er_select_glucosER (Solver& solver, unsigned int numClauses);
#endif

    ///// [ user_er_add_ ] /////
    // Description:
    //   User functions for introducing extension variables. These are heuristics for defining good extension variables.
    //
    // Parameters:
    //   solver          : the solver to add clauses to
    //   candidateClauses: clauses to use when selecting extended variable definitions
    //   maxNumNewVars   : the maximum number of new extension variables to introduce
    //
    // Return value:
    //   The function should return a map with the new extension variable as the key, and the pair of literals it corresponds
    //   to as the value. For simplicity, we restrict the interface to use Tseitin's original version of extension variables.
    //   Extension variables are defined as equivalent to a disjunction of a pair of literals. Complex extension variable
    //   definitions must be encoded with multiple extension variables. The size of the map should equal the number of new
    //   extension variables.

    static std::vector< std::pair<Lit, Lit> > getFreqSubexprs(std::tr1::unordered_map<std::pair<Lit, Lit>, int>& subexpr_counts, Solver& solver, unsigned int numSubexprs);

#if ER_USER_ADD_HEURISTIC != ER_ADD_HEURISTIC_NONE
    static std::vector< std::pair< Var, std::pair<Lit, Lit> > > user_er_add(Solver& solver, std::vector<CRef>& candidateClauses, unsigned int maxNumNewVars);
#endif

#if ER_USER_ADD_HEURISTIC == ER_ADD_HEURISTIC_SUBEXPR
    // Subexpression-based literal selection - select the disjunction of literals which occurs the most often together.
    static std::vector< std::pair< Var, std::pair<Lit, Lit> > > user_er_add_subexpr(Solver& solver, std::vector<CRef>& candidateClauses, unsigned int maxNumNewVars);

#elif ER_USER_ADD_HEURISTIC == ER_ADD_HEURISTIC_RANDOM
    // Random literal selection - select two literals at random and define a new extension variable over them.
    static std::vector< std::pair< Var, std::pair<Lit, Lit> > > user_er_add_random(Solver& solver, std::vector<CRef>& candidateClauses, unsigned int maxNumNewVars);

#elif ER_USER_ADD_HEURISTIC == ER_ADD_HEURISTIC_GLUCOSER
    // GlucosER literal selection - select the two literals which are not shared and define a new extension variable over them.
    static std::vector< std::pair< Var, std::pair<Lit, Lit> > > user_er_add_glucosER(Solver& solver, std::vector<CRef>& candidateClauses, unsigned int maxNumNewVars);
#endif

    ///// [ user_er_delete_ ] /////
    // Description:
    //   User functions for deleting extension variables. These are heuristics for deleting bad extension variables.
    //
    // Parameters:
    //   solver : the solver to remove extension variables from
    //
    // Return value:
    //   The function should return a list of extension variables that should be deleted

    // Delete all extension variables
    static std::tr1::unordered_set<Var> user_er_delete_all(Solver& solver);

    // Delete extension variables with activity below a threshold
    static std::tr1::unordered_set<Var> user_er_delete_activity(Solver& solver);

    // EXTENDED RESOLUTION - statistics
    // Functions for measuring extended resolution overhead
    void   extTimerStart();
    void   extTimerStop(struct rusage& ext_overhead);

    // Static helpers:
    //

    // Returns a random float 0 <= x < 1. Seed must never be 0.
    static inline double drand(double& seed) {
        seed *= 1389796;
        int q = (int)(seed / 2147483647);
        seed -= (double)q * 2147483647;
        return seed / 2147483647; }

    // Returns a random integer 0 <= x < size. Seed must never be 0.
    static inline int irand(double& seed, int size) {
        return (int)(drand(seed) * size); }
};


//=================================================================================================
// Implementation of inline methods:

inline CRef Solver::reason(Var x) const { return vardata[x].reason; }
inline int  Solver::level (Var x) const { return vardata[x].level; }

inline void Solver::insertVarOrder(Var x) {
    if (!order_heap.inHeap(x) && decision[x]) order_heap.insert(x); }

#if BRANCHING_HEURISTIC == VSIDS
inline void Solver::varDecayActivity() { var_inc *= (1 / var_decay); }
inline void Solver::varBumpActivity(Var v) { varBumpActivity(v, var_inc); }
inline void Solver::varBumpActivity(Var v, double inc) {
    if ( (activity[v] += inc) > 1e100 ) {
        // Rescale:
        for (int i = 0; i < nVars(); i++)
            activity[i] *= 1e-100;
        var_inc *= 1e-100; }

    // Update order_heap with respect to new activity:
    if (order_heap.inHeap(v))
        order_heap.decrease(v); }
#endif
#if ! LBD_BASED_CLAUSE_DELETION
inline void Solver::claDecayActivity() { cla_inc *= (1 / clause_decay); }
inline void Solver::claBumpActivity (Clause& c) {
        if ( (c.activity() += cla_inc) > 1e20 ) {
            // Rescale:
            for (int i = 0; i < learnts.size(); i++)
                ca[learnts[i]].activity() *= 1e-20;
            cla_inc *= 1e-20; } }
#endif

inline void Solver::checkGarbage(void){ return checkGarbage(garbage_frac); }
inline void Solver::checkGarbage(double gf){
    if (ca.wasted() > ca.size() * gf)
        garbageCollect(); }

// NOTE: enqueue does not set the ok flag! (only public methods do)
inline bool     Solver::enqueue         (Lit p, CRef from)                    { return value(p) != l_Undef ? value(p) != l_False : (uncheckedEnqueue(p, from), true); }
inline bool     Solver::addClause_      (      vec<Lit>& ps)                  { return addClauseToDB(clauses, ps); }
inline bool     Solver::addClause       (const vec<Lit>& ps)                  { ps.copyTo(add_tmp); return addClause_(add_tmp); }
inline bool     Solver::addEmptyClause  ()                                    { add_tmp.clear(); return addClause_(add_tmp); }
inline bool     Solver::addClause       (Lit p)                               { add_tmp.clear(); add_tmp.push(p); return addClause_(add_tmp); }
inline bool     Solver::addClause       (Lit p, Lit q)                        { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); return addClause_(add_tmp); }
inline bool     Solver::addClause       (Lit p, Lit q, Lit r)                 { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); add_tmp.push(r); return addClause_(add_tmp); }
inline bool     Solver::addClauseToDB   (vec<CRef>& db, Lit p)                { add_tmp.clear(); add_tmp.push(p); return addClauseToDB(db, add_tmp);}
inline bool     Solver::addClauseToDB   (vec<CRef>& db, Lit p, Lit q)         { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); return addClauseToDB(db, add_tmp);}
inline bool     Solver::addClauseToDB   (vec<CRef>& db, Lit p, Lit q, Lit r)  { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); add_tmp.push(r); return addClauseToDB(db, add_tmp);}
inline bool     Solver::locked          (const Clause& c) const               { return value(c[0]) == l_True && reason(var(c[0])) != CRef_Undef && ca.lea(reason(var(c[0]))) == &c; }
inline void     Solver::newDecisionLevel()                                    { trail_lim.push(trail.size()); }

// EXTENDED RESOLUTION
inline int      Solver::isExtVar     (Var x)           const { return (x >= originalNumVars) ? 1 : 0; }
inline int      Solver::getNumExtVars(const Clause& c) const {
    int numExtVar = 0;
    for (int i = 0; i < c.size(); i++) numExtVar += isExtVar(var(c[i]));
    return numExtVar;
}

inline int      Solver::decisionLevel ()      const   { return trail_lim.size(); }
inline uint32_t Solver::abstractLevel (Var x) const   { return 1 << (level(x) & 31); }
inline lbool    Solver::value         (Var x) const   { return assigns[x]; }
inline lbool    Solver::value         (Lit p) const   { return assigns[var(p)] ^ sign(p); }
inline lbool    Solver::modelValue    (Var x) const   { return model[x]; }
inline lbool    Solver::modelValue    (Lit p) const   { return model[var(p)] ^ sign(p); }
inline int      Solver::nAssigns      ()      const   { return trail.size(); }
inline int      Solver::nClauses      ()      const   { return clauses.size(); }
inline int      Solver::nLearnts      ()      const   { return learnts.size(); }
inline int      Solver::nVars         ()      const   { return vardata.size(); }

// EXTENDED RESOLUTION - statistics
inline int      Solver::nExtLearnts   ()      const   { return extLearnts.size(); }
inline int      Solver::nExtDefs      ()      const   { return extDefs.size(); }
inline int      Solver::nExtVars      ()      const   { return vardata.size() - originalNumVars; }

inline int      Solver::nFreeVars     ()      const   { return (int)dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]); }
inline void     Solver::setPolarity   (Var v, bool b) { polarity[v] = b; }
inline void     Solver::setDecisionVar(Var v, bool b) 
{ 
    if      ( b && !decision[v]) dec_vars++;
    else if (!b &&  decision[v]) dec_vars--;

    decision[v] = b;
    insertVarOrder(v);
}
inline void     Solver::setConfBudget(int64_t x){ conflict_budget    = conflicts    + x; }
inline void     Solver::setPropBudget(int64_t x){ propagation_budget = propagations + x; }
inline void     Solver::interrupt(){ asynch_interrupt = true; }
inline void     Solver::clearInterrupt(){ asynch_interrupt = false; }
inline void     Solver::budgetOff(){ conflict_budget = propagation_budget = -1; }
inline bool     Solver::withinBudget() const {
    return !asynch_interrupt &&
           (conflict_budget    < 0 || conflicts < (uint64_t)conflict_budget) &&
           (propagation_budget < 0 || propagations < (uint64_t)propagation_budget); }

inline bool     Solver::interrupted() const { return asynch_interrupt; }

// FIXME: after the introduction of asynchronous interrruptions the solve-versions that return a
// pure bool do not give a safe interface. Either interrupts must be possible to turn off here, or
// all calls to solve must return an 'lbool'. I'm not yet sure which I prefer.
inline bool     Solver::solve         ()                    { budgetOff(); assumptions.clear(); return solve_() == l_True; }
inline bool     Solver::solve         (Lit p)               { budgetOff(); assumptions.clear(); assumptions.push(p); return solve_() == l_True; }
inline bool     Solver::solve         (Lit p, Lit q)        { budgetOff(); assumptions.clear(); assumptions.push(p); assumptions.push(q); return solve_() == l_True; }
inline bool     Solver::solve         (Lit p, Lit q, Lit r) { budgetOff(); assumptions.clear(); assumptions.push(p); assumptions.push(q); assumptions.push(r); return solve_() == l_True; }
inline bool     Solver::solve         (const vec<Lit>& assumps){ budgetOff(); assumps.copyTo(assumptions); return solve_() == l_True; }
inline lbool    Solver::solveLimited  (const vec<Lit>& assumps){ assumps.copyTo(assumptions); return solve_(); }
inline bool     Solver::okay          ()      const   { return ok; }

inline void     Solver::toDimacs     (const char* file){ vec<Lit> as; toDimacs(file, as); }
inline void     Solver::toDimacs     (const char* file, Lit p){ vec<Lit> as; as.push(p); toDimacs(file, as); }
inline void     Solver::toDimacs     (const char* file, Lit p, Lit q){ vec<Lit> as; as.push(p); as.push(q); toDimacs(file, as); }
inline void     Solver::toDimacs     (const char* file, Lit p, Lit q, Lit r){ vec<Lit> as; as.push(p); as.push(q); as.push(r); toDimacs(file, as); }


//=================================================================================================
// Debug etc:

// EXTENDED RESOLUTION - statistics
inline void Solver::extTimerStart() {
    getrusage(RUSAGE_SELF, &ext_timer_start);
}
inline void Solver::extTimerStop(struct rusage& ext_overhead) {
    getrusage(RUSAGE_SELF, &ext_timer_end);
    
    // Add to total overhead
    ext_overhead.ru_utime.tv_sec  += ext_timer_end.ru_utime.tv_sec - ext_timer_start.ru_utime.tv_sec;
    ext_overhead.ru_utime.tv_usec += ext_timer_end.ru_utime.tv_usec;

    // Check if subtracting the initial time would result in underflow
    if (ext_timer_start.ru_utime.tv_usec > ext_overhead.ru_utime.tv_usec) {
        ext_overhead.ru_utime.tv_usec += 1000000;
        ext_overhead.ru_utime.tv_sec  -= 1;
    }
    ext_overhead.ru_utime.tv_usec -= ext_timer_start.ru_utime.tv_usec;

    // Check if we carry over to the next second
    if (ext_overhead.ru_utime.tv_usec >= 1000000) {
        ext_overhead.ru_utime.tv_usec -= 1000000;
        ext_overhead.ru_utime.tv_sec  += 1;
    }
}

static inline double readTimer(struct rusage& ext_overhead) {
    return (double)ext_overhead.ru_utime.tv_sec + (double)ext_overhead.ru_utime.tv_usec / 1000000;
}
inline double Solver::extTimerRead(unsigned int i) {
    switch(i) {
        case 0: return readTimer(ext_sel_overhead);
        case 1: return readTimer(ext_add_overhead);
        case 2: return readTimer(ext_delC_overhead);
        case 3: return readTimer(ext_delV_overhead);
        case 4: return readTimer(ext_sub_overhead);
        case 5: return readTimer(ext_stat_overhead);
        default: return -1.;
    }
}

//=================================================================================================
}

#endif
