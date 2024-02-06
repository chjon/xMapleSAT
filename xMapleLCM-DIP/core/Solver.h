/****************************************************************************************[Solver.h]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010, Niklas Sorensson

Chanseok Oh's MiniSat Patch Series -- Copyright (c) 2015, Chanseok Oh
 
Maple_LCM, Based on MapleCOMSPS_DRUP -- Copyright (c) 2017, Mao Luo, Chu-Min LI, Fan Xiao: implementing a learnt clause minimisation approach
Reference: M. Luo, C.-M. Li, F. Xiao, F. Manya, and Z. L. , “An effective learnt clause minimization approach for cdcl sat solvers,” in IJCAI-2017, 2017, pp. to–appear.
 
Maple_LCM_Dist, Based on Maple_LCM -- Copyright (c) 2017, Fan Xiao, Chu-Min LI, Mao Luo: using a new branching heuristic called Distance at the beginning of search
 
 
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

#define GLUCOSE23

#ifdef GLUCOSE23
#define INT_QUEUE_AVG
#define LOOSE_PROP_STAT
#endif

#define BRANCHING_HEURISTIC_DYNAMIC  0
#define BRANCHING_HEURISTIC_VSIDS    1
#define BRANCHING_HEURISTIC_CHB      2
#define BRANCHING_HEURISTIC_DISTANCE 3

#ifndef BRANCHING_HEURISTIC
#define BRANCHING_HEURISTIC BRANCHING_HEURISTIC_DYNAMIC
#endif

#include "core/AssignmentTrail.h"
#include "core/BranchingHeuristicManager.h"
#include "core/ClauseDatabase.h"
#include "core/ConflictAnalyzer.h"
#include "core/ProofLogger.h"
#include "core/PropagationQueue.h"
#include "core/RandomNumberGenerator.h"
#include "core/RestartHeuristicManager.h"
#include "core/SolverTypes.h"
#include "core/UnitPropagator.h"
#include "mtl/Alg.h"
#include "mtl/Vec.h"
#include "utils/Options.h"

namespace Minisat {

/**
 * @brief This is the main SAT solver class. It is responsible for solving instances of the SAT
 * problem.
 * 
 */
class Solver {
protected:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // MEMBER VARIABLES

    /// @brief If FALSE, the constraints are already unsatisfiable. No part of the solver state may
    /// be used!
    bool ok;

    /// @brief A flag that is set to true by the interrupt handler if the solver receives an exit
    /// signal. Used to clean up and exit gracefully.
    bool asynch_interrupt;

protected:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TEMPORARY VARIABLES

    /// @brief Used to reduce reallocation overhead associated with learning a new clause
    vec<Lit> learnt_clause;

    // Used to reduce reallocation overhead associated with proof generation
    vec<Lit> add_oc;

    /// @brief Used to reduce reallocation overhead associated with simplifying a clause
    vec<Lit> simp_learnt_clause;

    /// @brief Used to reduce reallocation overhead associated with simplifying a clause
    vec<CRef> simp_reason_clause;

protected:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // PARAMETERS

    /// @brief Current set of assumptions provided to solve by the user.
    vec<Lit> assumptions;

    /// @brief Resource constraint: the number of conflicts to search for before exiting. -1 means
    /// no budget.
    int64_t conflict_budget;

    // adjust simplifyAll occasion

    // Clause simplification frequency parameters

    /// @brief the threshold conflict fraction for triggering mass simplification upon restart
    long curSimplify;

    /// @brief Number of conflicts before next call to simplifyAll
    int nbconfbeforesimplify;

    /// @brief Increment for @code{nbconfbeforesimplify}
    int incSimplify;

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // PUBLIC API

    /// @brief The verbosity level of the solver. (0=silent, 1=some, 2=more)
    int verbosity;

    /// @brief Says whether the solver will produce a DRAT proof in file proof.txt
    int produce_proof;
    
    /// @brief If problem is satisfiable, this vector contains the model (if any).
    vec<lbool> model;

    /// @brief if problem is unsatisfiable (possibly under assumptions), this vector represents the
    /// final conflict clause expressed in the assumptions.
    vec<Lit> conflict;

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // STATISTICS

    /// @brief The number of calls to @code{solve_()}
    uint64_t solves;
    
    /// @brief The number of calls to @code{search()}
    uint64_t starts;
    
    /// @brief The total number of conflicts encountered during search.
    uint64_t conflicts;

    /// @brief The total number of conflicts with a DIP found encountered during search.
    uint64_t dip_conflicts;

    /// @brief Number of top-level assignments since last execution of 'simplify()'.
    int simpDB_assigns;

    /// @brief Remaining number of propagations that must be made before the next execution of
    /// 'simplify()'.
    int64_t simpDB_props;

    
public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // SOLVER COMPONENTS

    /// @brief Memory manager for allocating/deallocating clauses
    ClauseAllocator ca;

    /// @brief Assignment trail
    AssignmentTrail assignmentTrail;

    /// @brief Branching heuristic
    BranchingHeuristicManager branchingHeuristicManager;

    /// @brief Clause database
    ClauseDatabase clauseDatabase;

    /// @brief Conflict analyzer
    ConflictAnalyzer conflictAnalyzer;

    /// @brief Proof generator
    ProofLogger proofLogger;

    /// @brief Unit propagation queue
    PropagationQueue propagationQueue;

    /// @brief Random number generator
    RandomNumberGenerator randomNumberGenerator;

    /// @brief Restart heuristic
    RestartHeuristicManager restartHeuristicManager;

    /// @brief Unit propagator
    UnitPropagator unitPropagator;

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // CONSTRUCTORS

    /**
     * @brief Construct a new Solver object
     * 
     */
    Solver();

    /**
     * @brief Destroy the Solver object
     * 
     */
    virtual ~Solver() = default;

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // ACCESSORS

    /**
     * @brief Check whether the solver is in a conflicting state
     * 
     * @return false iff the solver is in a conflicting state
     */
    bool okay(void) const;

    /**
     * @brief The value of a variable in the last model.
     * 
     * @param x the variable whose value should be returned.
     * @return the truth assignment of the variable in the last model.
     * 
     * @pre The last call to solve must have been satisfiable.
     */
    lbool modelValue(Var x) const;

    /**
     * @brief The value of a literal in the last model.
     * 
     * @param x the literal whose value should be returned.
     * @return the truth assignment of the literal in the last model.
     * 
     * @pre The last call to solve must have been satisfiable.
     */
    lbool modelValue(Lit p) const;

    /**
     * @brief Check whether the solver has been interrupted
     * 
     * @return true iff @code{async_interrupt} is true
     */
    bool interrupted(void) const;

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // PARAMETER MODIFICATION

    /**
     * @brief Set the conflict budget
     * 
     * @param x the remaining number of conflicts.
     */
    void setConfBudget(int64_t x);

    /**
     * @brief Turn off the conflict budget.
     * 
     */
    void budgetOff(void);

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // STATE MODIFICATION

    /**
     * @brief Relocate all clauses
     * 
     * @param to the ClauseAllocator to relocate to
     */
    virtual void relocAll(ClauseAllocator& to);

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // PROBLEM SPECIFICATION

    /**
     * @brief Creates a new SAT variable in the solver
     * 
     * @param polarity The preferred initial polarity of the new variable
     * @param dvar true iff the new variable should be used as a decision variable (NOTE! This has
     * effects on the meaning of a SATISFIABLE result).
     * @return the ID of the new variable
     */
    virtual Var newVar (bool polarity = true, bool dvar = true);
    
    /**
     * @brief Add a new input clause
     * 
     * @param ps the literals of the new clause
     * @return true iff the solver is in a consistent state after adding the clause
     */
    bool addClause(vec<Lit>& ps);

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // SOLVING

    /**
     * @brief Simplify the clause database according to the current top-level assigment.
     * 
     * @details Currently, the only thing done here is the removal of satisfied clauses, but more
     * things can be put here.
     * 
     * @return false if the solver is already in a conflicting state, true otherwise
     */
    virtual bool simplify(void);

    bool simplifyAll();

    /**
     * @brief Search for a model that respects a given set of assumptions.
     * 
     * @param assumps a set of assumptions
     * @return true if the formula is found to be satisfiable, false otherwise
     */
    bool solve(const vec<Lit>& assumps);

    /**
     * @brief Search for a model that respects a given set of assumptions (With resource constraints).
     * 
     * @param assumps a set of assumptions
     * @return l_True if the formula is satisfiable, l_False if the formula is unsatisfiable, or
     * l_Undef if the satisfiability of the formula is unknown.
     */
    lbool solveLimited(const vec<Lit>& assumps);

    /**
     * @brief Search without assumptions.
     * 
     * @return true if the formula is found to be satisfiable, false otherwise
     */
    bool solve(void);

    /**
     * @brief Search for a model that respects a single assumption.
     * 
     * @param p the first assumption
     * @return true if the formula is found to be satisfiable, false otherwise
     */
    bool solve(Lit p);

    /**
     * @brief Search for a model that respects two assumptions.
     * 
     * @param p the first assumption
     * @param q the second assumption
     * @return true if the formula is found to be satisfiable, false otherwise 
     */
    bool solve(Lit p, Lit q);

    /**
     * @brief Search for a model that respects three assumptions.
     * 
     * @param p the first assumption
     * @param q the second assumption
     * @param r the third assumption
     * @return true if the formula is found to be satisfiable, false otherwise 
     */
    bool solve(Lit p, Lit q, Lit r);

public:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // EVENT HANDLERS

    /**
     * @brief Trigger a (potentially asynchronous) interruption of the solver.
     * 
     */
    void interrupt(void);

    /**
     * @brief Clear interrupt indicator flag.
     * 
     */
    void clearInterrupt(void);

    /**
     * @brief Update data structures when a clause is deleted
     * 
     * @param c the clause that was deleted
     * @param cr the reference to the clause that was deleted
     */
    virtual void handleEventClauseDeleted(const Clause& c, CRef cr);

protected:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // HELPER FUNCTIONS

    /**
     * @brief Search for a given number of conflicts
     * 
     * @param nof_conflicts The number of conflicts to search for before exiting the function
     * @return 'l_True' if a partial assigment that is consistent with respect to the clause set
     * is found. If all variables are decision variables, this means that the clause set is
     * satisfiable. 'l_False' if the clause set is unsatisfiable. 'l_Undef' if the bound on number
     * of conflicts is reached.
     */
    virtual lbool search(int& nof_conflicts);
    
    /**
     * @brief Main solve method (assumptions given in 'assumptions')
     * 
     * @return l_True if the formula is satisfiable, l_False if the formula is unsatisfiable, or
     * l_Undef if the satisfiability of the formula is unknown.
     */
    virtual lbool solve_(void);

    /**
     * @brief Check whether the solver should exit searching early
     * 
     * @return true if the solver should stop searching.
     */
    bool withinBudget(void) const;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // PRINT CLAUSE INFORMATION
  
    public:
      /**
       * @brief Print clause with literal values and decision levels
       * 
       * @param c the clause that is to be printed
       * @param out the stream to which the output is redirected
       */
      
      void writeClauseWithValues (const Clause& c, ostream& out); 
      void writeClauseWithValues (const vector<Lit>& c, ostream& out);
      void writeClauseWithValues (const vec<Lit>& c, ostream& out);

  
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // HELPER FUNCTIONS FOR simplifyAll

    template<int db_mark>
    bool simplifyLearntDB();

    void removeFalseLits(Clause& c);

    bool detachAndSimplify(CRef cr);

    void simplifyLearnt(Clause& c);

private:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // FRIENDS

    friend AssignmentTrail;
    friend ClauseDatabase;
    friend UnitPropagator;
    friend BranchingHeuristicManager;
    friend ConflictAnalyzer;
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION OF INLINE METHODS

//////////////
// ACCESSORS

inline bool Solver::okay(void) const {
    return ok;
}

inline lbool Solver::modelValue(Var x) const {
    return model[x];
}

inline lbool Solver::modelValue(Lit p) const {
    return model[var(p)] ^ sign(p);
}

inline bool Solver::interrupted(void) const {
    return asynch_interrupt;
}

///////////////////////////
// PARAMETER MODIFICATION

inline void Solver::setConfBudget(int64_t x) {
    conflict_budget = conflicts + x;
}

inline void Solver::budgetOff() {
    conflict_budget = -1;
    unitPropagator.budgetOff();
}

///////////////////////
// STATE MODIFICATION

inline void Solver::relocAll(ClauseAllocator& to) {
    unitPropagator .relocAll(to);
    assignmentTrail.relocAll(to);
    clauseDatabase .relocAll(to);
}

//////////////////////////
// PROBLEM SPECIFICATION

inline Var Solver::newVar(bool sign, bool dvar) {
    int v = assignmentTrail  .newVar();
    clauseDatabase           .newVar(v);
    propagationQueue         .newVar(v);
    unitPropagator           .newVar(v);
    branchingHeuristicManager.newVar(v, sign, dvar);
    conflictAnalyzer         .newVar(v);
    return v;
}

///////////////////
// EVENT HANDLERS

inline void Solver::interrupt() {
    asynch_interrupt = true;
}

inline void Solver::clearInterrupt() {
    asynch_interrupt = false;
}

inline void Solver::handleEventClauseDeleted(const Clause& c, CRef cr) {
    // Log clause removal
    proofLogger.removeClause(c);

    // Don't leave pointers to free'd memory!
    unitPropagator.detachClause(cr, false);
    assignmentTrail.handleEventClauseDeleted(c);
}

////////////
// SOLVING

// FIXME: after the introduction of asynchronous interrruptions the solve-versions that return a
// pure bool do not give a safe interface. Either interrupts must be possible to turn off here, or
// all calls to solve must return an 'lbool'. I'm not yet sure which I prefer.
inline bool     Solver::solve         ()                    { budgetOff(); assumptions.clear(); return solve_() == l_True; }
inline bool     Solver::solve         (Lit p)               { budgetOff(); assumptions.clear(); assumptions.push(p); return solve_() == l_True; }
inline bool     Solver::solve         (Lit p, Lit q)        { budgetOff(); assumptions.clear(); assumptions.push(p); assumptions.push(q); return solve_() == l_True; }
inline bool     Solver::solve         (Lit p, Lit q, Lit r) { budgetOff(); assumptions.clear(); assumptions.push(p); assumptions.push(q); assumptions.push(r); return solve_() == l_True; }
inline bool     Solver::solve         (const vec<Lit>& assumps){ budgetOff(); assumps.copyTo(assumptions); return solve_() == l_True; }
inline lbool    Solver::solveLimited  (const vec<Lit>& assumps){ assumps.copyTo(assumptions); return solve_(); }

/////////////////////
// HELPER FUNCTIONS

inline bool Solver::withinBudget() const {
    return !asynch_interrupt &&
        (conflict_budget    < 0 || conflicts < (uint64_t)conflict_budget) &&
        unitPropagator.withinBudget();
}

/////////////////////////////////////
// HELPER FUNCTIONS FOR simplifyAll

inline void Solver::removeFalseLits(Clause& c) {
    int i, j;
    for (i = j = 0; i < c.size(); i++) {
        if (assignmentTrail.value(c[i]) != l_False)
            c[j++] = c[i];
    }
    c.shrink(i - j);
}

inline bool Solver::detachAndSimplify(CRef cr) {
  Clause& c = ca[cr];

    // Remove satisfied clauses and check for false literals
    bool false_lit = false;
    for (int i = 0; i < c.size(); i++){
        if (assignmentTrail.value(c[i]) == l_True){
            clauseDatabase.removeClause(cr);
            return true;
        } else if (assignmentTrail.value(c[i]) == l_False){
            false_lit = true;
        }
    }

    // Detach clause so it can be modified
    clauseDatabase.detachClause(cr, true);
    
    // Remove false literals
    int saved_size = c.size();
    if (false_lit) removeFalseLits(c);
    assert(c.size() > 1);

    // Further simplify the learnt clause
    simplifyLearnt(c);
    assert(c.size() > 0);
    
    // Log clause simplification
    if (saved_size != c.size()) proofLogger.addClause(c);

    return false;
}
}

#endif
