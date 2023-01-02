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
#include "core/PropagationQueue.h"
#include "core/UnitPropagator.h"
#include "core/BranchingHeuristicManager.h"
#include "core/ClauseDatabase.h"
#include "core/ConflictAnalyzer.h"

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
    Var newVar (bool polarity = true, bool dvar = true); // Add a new variable with parameters specifying variable mode.
    
    /**
     * @brief Add a new input clause
     * 
     * @param ps the literals of the new clause
     * @return true iff the solver is in a consistent state after adding the clause
     */
    bool addClause(vec<Lit>& ps);

    // Solving:
    //
    bool  simplify    ();                        // Removes already satisfied clauses.
    bool  solve       (const vec<Lit>& assumps); // Search for a model that respects a given set of assumptions.
    lbool solveLimited(const vec<Lit>& assumps); // Search for a model that respects a given set of assumptions (With resource constraints).
    bool  solve       ();                        // Search without assumptions.
    bool  solve       (Lit p);                   // Search for a model that respects a single assumption.
    bool  solve       (Lit p, Lit q);            // Search for a model that respects two assumptions.
    bool  solve       (Lit p, Lit q, Lit r);     // Search for a model that respects three assumptions.
    bool  okay        () const;                  // FALSE means solver is in a conflicting state

    // Read state:
    //
    lbool   modelValue (Var x) const; // The value of a variable in the last model. The last call to solve must have been satisfiable.
    lbool   modelValue (Lit p) const; // The value of a literal in the last model. The last call to solve must have been satisfiable.
    int     nFreeVars  ()      const;

    // Resource contraints:
    //
    void    setConfBudget(int64_t x);
    void    budgetOff();
    void    interrupt();          // Trigger a (potentially asynchronous) interruption of the solver.
    void    clearInterrupt();     // Clear interrupt indicator flag.

    // Extra results: (read-only member variable)
    //
    vec<lbool> model;    // If problem is satisfiable, this vector contains the model (if any).

    /// @brief if problem is unsatisfiable (possibly under assumptions), this vector represents the
    // final conflict clause expressed in the assumptions.
    vec<Lit> conflict;

    // Mode of operation:
    //
    int       verbosity;
#if ! LBD_BASED_CLAUSE_DELETION
    double    clause_decay;
#endif
    bool      luby_restart;

    int       restart_first;      // The initial restart limit.                                                                (default 100)
    double    restart_inc;        // The factor with which the restart limit is multiplied in each restart.                    (default 1.5)
    double    learntsize_factor;  // The intitial limit for learnt clauses is a factor of the original clauses.                (default 1 / 3)
    double    learntsize_inc;     // The limit for learnt clauses is multiplied with this factor each restart.                 (default 1.1)

    int       learntsize_adjust_start_confl;
    double    learntsize_adjust_inc;

    // Statistics: (read-only member variable)
    //
    uint64_t solves, starts, conflicts;

    uint64_t lbd_calls;
    vec<uint64_t> lbd_seen;

protected:

    // Solver state:
    //
    bool                ok;               // If FALSE, the constraints are already unsatisfiable. No part of the solver state may be used!
#if ! LBD_BASED_CLAUSE_DELETION
    double              cla_inc;          // Amount to bump next clause with.
#endif

    int                 simpDB_assigns;   // Number of top-level assignments since last execution of 'simplify()'.
    int64_t             simpDB_props;     // Remaining number of propagations that must be made before next execution of 'simplify()'.
    vec<Lit>            assumptions;      // Current set of assumptions provided to solve by the user.
    double              progress_estimate;// Set by 'search()'.

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

    // Maintaining Variable/Clause activity:
    //
#if ! LBD_BASED_CLAUSE_DELETION
    void     claDecayActivity ();                      // Decay all clauses with the specified factor. Implemented by increasing the 'bump' value instead.
    void     claBumpActivity  (Clause& c);             // Increase a clause with the current 'bump' value.
#endif

    // Misc:
    //
    bool withinBudget() const;

public:
    VariableDatabase          variableDatabase;
    ClauseAllocator           ca; // Memory manager for allocating/deallocating clauses
    RandomNumberGenerator     randomNumberGenerator;
    AssignmentTrail           assignmentTrail;
    PropagationQueue          propagationQueue;
    UnitPropagator            unitPropagator;
    BranchingHeuristicManager branchingHeuristicManager;
    ClauseDatabase            clauseDatabase;
    ConflictAnalyzer          conflictAnalyzer;

private:
    friend AssignmentTrail;
    friend ClauseDatabase;
    friend UnitPropagator;
    friend BranchingHeuristicManager;
    friend ConflictAnalyzer;
};


//=================================================================================================
// Implementation of inline methods:

#if ! LBD_BASED_CLAUSE_DELETION
inline void Solver::claDecayActivity() { cla_inc *= (1 / clause_decay); }
inline void Solver::claBumpActivity (Clause& c) {
        if ( (c.activity() += cla_inc) > 1e20 ) {
            // Rescale:
            for (int i = 0; i < learnts.size(); i++)
                ca[learnts[i]].activity() *= 1e-20;
            cla_inc *= 1e-20; } }
#endif

inline lbool    Solver::modelValue    (Var x) const   { return model[x]; }
inline lbool    Solver::modelValue    (Lit p) const   { return model[var(p)] ^ sign(p); }
inline int      Solver::nFreeVars     ()      const   { return (int)branchingHeuristicManager.dec_vars - assignmentTrail.nRootAssigns(); }
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

//=================================================================================================
// Debug etc:


//=================================================================================================
}

#endif
