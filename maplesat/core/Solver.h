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

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "mtl/Multiheap.h"
#include "mtl/Alg.h"
#include "utils/Options.h"
#include "core/SolverTypes.h"
#include "core/VariableDatabase.h"
#include "core/RandomNumberGenerator.h"
#include "core/AssignmentTrail.h"
#include "core/UnitPropagator.h"

namespace Minisat {

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
#if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
    lbool   bcpValue  (Var x) const;       // The queued value of a variable.
    lbool   bcpValue  (Lit p) const;       // The queued value of a literal.
#endif
    lbool   modelValue (Var x) const;       // The value of a variable in the last model. The last call to solve must have been satisfiable.
    lbool   modelValue (Lit p) const;       // The value of a literal in the last model. The last call to solve must have been satisfiable.
    int     nClauses   ()      const;       // The current number of original clauses.
    int     nLearnts   ()      const;       // The current number of learnt clauses.
    int     nFreeVars  ()      const;

    // Resource contraints:
    //
    void    setConfBudget(int64_t x);
    void    budgetOff();
    void    interrupt();          // Trigger a (potentially asynchronous) interruption of the solver.
    void    clearInterrupt();     // Clear interrupt indicator flag.

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
    bool      luby_restart;
    int       ccmin_mode;         // Controls conflict clause minimization (0=none, 1=basic, 2=deep).
    int       phase_saving;       // Controls the level of phase saving (0=none, 1=limited, 2=full).
    bool      rnd_pol;            // Use random polarities for branching heuristics.
    bool      rnd_init_act;       // Initialize variable activities with a small random value.
    double    garbage_frac;       // The fraction of wasted memory allowed before a garbage collection is triggered.

    int       restart_first;      // The initial restart limit.                                                                (default 100)
    double    restart_inc;        // The factor with which the restart limit is multiplied in each restart.                    (default 1.5)
    double    learntsize_factor;  // The intitial limit for learnt clauses is a factor of the original clauses.                (default 1 / 3)
    double    learntsize_inc;     // The limit for learnt clauses is multiplied with this factor each restart.                 (default 1.1)

    int       learntsize_adjust_start_confl;
    double    learntsize_adjust_inc;

    // Statistics: (read-only member variable)
    //
    uint64_t solves, starts, decisions, rnd_decisions, conflicts;
    uint64_t dec_vars, clauses_literals, learnts_literals, max_literals, tot_literals;

    uint64_t lbd_calls;
    vec<uint64_t> lbd_seen;
    vec<uint64_t> picked;
    vec<uint64_t> conflicted;
#if PRIORITIZE_ER
    // Number of times a variable occurs in a clause
    vec<uint64_t> degree;
    // Map from variables to their extension level
    vec<uint64_t> extensionLevel;

#ifdef EXTLVL_ACTIVITY
    vec<double> extensionLevelActivity; // The activity of each extension level
#endif
#ifdef POLARITY_VOTING
    vec<double> group_polarity;   // The preferred polarity of each group.
#endif
#endif
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

    struct VarOrderLt {
        const vec<double>&  activity;
#if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
        const vec<uint64_t>& order_param;
        bool greater_than;
        bool operator () (Var x, Var y) const {
            if (order_param[x] != order_param[y]) return greater_than ^ (order_param[x] < order_param[y]);
            return activity[x] > activity[y];
        }
        VarOrderLt(const vec<double>&  act, const vec<uint64_t>& ord, bool gt = false)
            : activity(act)
            , order_param(ord)
            , greater_than(gt)
        { }
#else
        bool operator () (Var x, Var y) const { return activity[x] > activity[y]; }
        VarOrderLt(const vec<double>&  act) : activity(act) { }
#endif
    };

    // Solver state:
    //
    bool                ok;               // If FALSE, the constraints are already unsatisfiable. No part of the solver state may be used!
    vec<CRef>           clauses;          // List of problem clauses.
    vec<CRef>           learnts;          // List of learnt clauses.
#if ! LBD_BASED_CLAUSE_DELETION
    double              cla_inc;          // Amount to bump next clause with.
#endif
    vec<double>         activity;         // A heuristic measurement of the activity of a variable.
    double              var_inc;          // Amount to bump next variable with.
    vec<char>           polarity;         // The preferred polarity of each variable.
    vec<char>           decision;         // Declares if a variable is eligible for selection in the decision heuristic.

    int                 simpDB_assigns;   // Number of top-level assignments since last execution of 'simplify()'.
    int64_t             simpDB_props;     // Remaining number of propagations that must be made before next execution of 'simplify()'.
    vec<Lit>            assumptions;      // Current set of assumptions provided to solve by the user.
#if PRIORITIZE_ER
#ifdef EXTLVL_ACTIVITY
    Multiheap<VarOrderLt> order_heap;
#else
    Heap<VarOrderLt>    order_heap_extlvl;       // A priority queue of variables ordered with respect to the extension level.
    Heap<VarOrderLt>    order_heap_degree;       // A priority queue of variables ordered with respect to the variable degree.
#endif
#else
    Heap<VarOrderLt>    order_heap;       // A priority queue of variables ordered with respect to the variable activity.
#endif
    double              progress_estimate;// Set by 'search()'.
    bool                remove_satisfied; // Indicates whether possibly inefficient linear scan for satisfied clauses should be performed in 'simplify'.

    ClauseAllocator     ca;

    // Temporaries (to reduce allocation overhead). Each variable is prefixed by the method in which it is
    // used, exept 'seen' wich is used in several places.
    //
    vec<char>           seen;
    vec<Lit>            analyze_stack;
    vec<Lit>            analyze_toclear;
    vec<Lit>            add_tmp;

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
    void     analyze          (CRef confl, vec<Lit>& out_learnt, int& out_btlevel);    // (bt = backtrack)
    void     analyzeFinal     (Lit p, vec<Lit>& out_conflict);                         // COULD THIS BE IMPLEMENTED BY THE ORDINARIY "analyze" BY SOME REASONABLE GENERALIZATION?
    bool     litRedundant     (Lit p, uint32_t abstract_levels);                       // (helper method for 'analyze()')

    template<class V> int lbd (const V& clause) {
        lbd_calls++;
        int lbd = 0;
        for (int i = 0; i < clause.size(); i++) {
            int l = assignmentTrail.level(var(clause[i]));
            if (lbd_seen[l] != lbd_calls) {
                lbd++;
                lbd_seen[l] = lbd_calls;
            }
        }
        return lbd;
    }
    lbool    search           (int nof_conflicts); // Search for a given number of conflicts.
    lbool    solve_           ();                  // Main solve method (assumptions given in 'assumptions').
    void     reduceDB         ();                  // Reduce the set of learnt clauses.
    void     removeSatisfied  (vec<CRef>& cs);     // Shrink 'cs' to contain only non-satisfied clauses.
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
    bool     satisfied        (const Clause& c) const; // Returns TRUE if a clause is satisfied in the current state.

    void     relocAll         (ClauseAllocator& to);

    // Misc:
    //
    bool     withinBudget     ()      const;

public:
    VariableDatabase      variableDatabase;
    RandomNumberGenerator randomNumberGenerator;
    AssignmentTrail       assignmentTrail;
    UnitPropagator        unitPropagator;

private:
    friend AssignmentTrail;
    friend UnitPropagator;
};


//=================================================================================================
// Implementation of inline methods:

inline void Solver::insertVarOrder(Var x) {
#if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
    Heap<VarOrderLt>& order_heap = extensionLevel[x] ? order_heap_extlvl : order_heap_degree;
#endif
    if (!order_heap.inHeap(x) && decision[x]) order_heap.insert(x); }

#if BRANCHING_HEURISTIC == VSIDS
inline void Solver::varDecayActivity() { var_inc *= (1 / var_decay); }
inline void Solver::varBumpActivity(Var v) { varBumpActivity(v, var_inc); }
inline void Solver::varBumpActivity(Var v, double inc) {
    if ( (activity[v] += inc) > 1e100 ) {
#if PRIORITIZE_ER && defined(EXTLVL_ACTIVITY)
        // Clear extension level activity
        for (int i = 0; i < extensionLevelActivity.size(); i++) {
            extensionLevelActivity[i] = 0;
        }
#endif
        // Rescale:
        for (int i = 0; i < variableDatabase.nVars(); i++) {
            activity[i] *= 1e-100;
#if PRIORITIZE_ER && defined(EXTLVL_ACTIVITY)
            extensionLevelActivity[extensionLevel[i]] += activity[i];
#endif
        }
        var_inc *= 1e-100;
    }

    // Update order_heap with respect to new activity:
#if PRIORITIZE_ER && !defined(EXTLVL_ACTIVITY)
    if (order_heap_extlvl.inHeap(v))
        order_heap_extlvl.decrease(v);
    if (order_heap_degree.inHeap(v))
        order_heap_degree.decrease(v);
#else
    if (order_heap.inHeap(v))
        order_heap.decrease(v);
#endif
}
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
inline bool     Solver::addClause       (const vec<Lit>& ps)    { ps.copyTo(add_tmp); return addClause_(add_tmp); }
inline bool     Solver::addEmptyClause  ()                      { add_tmp.clear(); return addClause_(add_tmp); }
inline bool     Solver::addClause       (Lit p)                 { add_tmp.clear(); add_tmp.push(p); return addClause_(add_tmp); }
inline bool     Solver::addClause       (Lit p, Lit q)          { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); return addClause_(add_tmp); }
inline bool     Solver::addClause       (Lit p, Lit q, Lit r)   { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); add_tmp.push(r); return addClause_(add_tmp); }

#if BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
inline lbool    Solver::bcpValue      (Var x) const   { return bcp_assigns[x]; }
inline lbool    Solver::bcpValue      (Lit p) const   { return bcp_assigns[var(p)] ^ sign(p); }
#endif
inline lbool    Solver::modelValue    (Var x) const   { return model[x]; }
inline lbool    Solver::modelValue    (Lit p) const   { return model[var(p)] ^ sign(p); }
inline int      Solver::nClauses      ()      const   { return clauses.size(); }
inline int      Solver::nLearnts      ()      const   { return learnts.size(); }
inline int      Solver::nFreeVars     ()      const   { return (int)dec_vars - assignmentTrail.nRootAssigns(); }
inline void     Solver::setPolarity   (Var v, bool b) { polarity[v] = b; }
inline void     Solver::setDecisionVar(Var v, bool b) { 
    if      ( b && !decision[v]) dec_vars++;
    else if (!b &&  decision[v]) dec_vars--;

    decision[v] = b;
    insertVarOrder(v);
}
inline void     Solver::setConfBudget(int64_t x){ conflict_budget    = conflicts    + x; }
inline void     Solver::interrupt(){ asynch_interrupt = true; }
inline void     Solver::clearInterrupt(){ asynch_interrupt = false; }
inline void     Solver::budgetOff(){ conflict_budget = propagation_budget = -1; }
inline bool     Solver::withinBudget() const {
    return !asynch_interrupt &&
        (conflict_budget    < 0 || conflicts < (uint64_t)conflict_budget) &&
        unitPropagator.withinBudget();
}

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


//=================================================================================================
}

#endif
