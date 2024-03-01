/******************************************************************************[ConflictAnalyzer.h]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010, Niklas Sorensson

MapleSAT_Refactor, based on MapleSAT -- Copyright (c) 2022, Jonathan Chung, Vijay Ganesh, Sam Buss

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

#ifndef Minisat_ConflictAnalyzer_h
#define Minisat_ConflictAnalyzer_h

#include "core/SolverTypes.h"
#include "core/TwoVertexBottlenecks.h"

#include "er/ERManager.h"

#include <iostream>
#include <stack>
using namespace std;

namespace Minisat {
  // Forward declarations
  class Solver;
  class AssignmentTrail;
  class BranchingHeuristicManager;

  // DIPGraphEncoder maps solver literals to integers used in the DIP-detection algorithm
  class DIPGraphEncoder {

    int DIP_ENCODER_UNDEF = -1;
    int nextVertex; // vertices start from 0 (conflict node)
    vec<int> _solver2sam;
    vec<Lit> _sam2solver;

  public:
    int numVertices ( ) const;
    DIPGraphEncoder (int nVars ):
      nextVertex(0), _solver2sam(nVars+1,DIP_ENCODER_UNDEF) { }
    Lit Sam2Solver (int var) const;
    int Solver2Sam (Lit lit);
  };
  
  inline int DIPGraphEncoder::numVertices ( )      const { return nextVertex; }
  inline Lit DIPGraphEncoder::Sam2Solver (int var) const { return _sam2solver[var]; };
  inline int DIPGraphEncoder::Solver2Sam (Lit lit)       {
    int v = var(lit);
    if (_solver2sam[v] == DIP_ENCODER_UNDEF) {
      _solver2sam[v] = nextVertex;
      _sam2solver.push(lit);
      ++nextVertex;
      //cout << "Lit " << lit << " gets " << nextVertex - 1 << std::endl;
      return nextVertex-1;
    }
    else return _solver2sam[v];
  };

  /**
   * @brief This class is responsible for analyzing conflict graphs to generate and simplify
   * learnt clauses.
   * 
   */
  class ConflictAnalyzer {
  protected:
    ///////////////////////////////////////////////////////////////////////////////////////////
    // HELPER TYPES

    enum ConflictClauseMinimizationMode: int {
      NONE  = 0,
      BASIC = 1,
      DEEP  = 2,
    };

  private:
    ///////////////////////////////////////////////////////////////////////////////////////////
    // SOLVER REFERENCES

    AssignmentTrail& assignmentTrail;
    BranchingHeuristicManager& branchingHeuristicManager;
    ClauseAllocator& ca;
    Solver& solver;
  protected:
    ///////////////////////////////////////////////////////////////////////////////////////////
    // PARAMETERS
        
    /// @brief Controls conflict clause minimization
    ConflictClauseMinimizationMode ccmin_mode;

    // @brief Controls computation of DIP
    bool compute_dip;

    // @brief Minimum number of times a DIP has to appear before we introduce it
    int dip_pair_threshold;

    // @brief Window size to determine whether DIP is better than average over sliding window
    int dip_window_size;

    // @brief Introduces DIP only if clause DIP -> conflict has LBD <= dip_max_LBD
    int dip_max_LBD;

    static const int MIDDLE_DIP = 1;
    static const int CLOSEST_TO_CONFLICT = 2;
    static const int RANDOM_DIP = 3;

    // @brief Type of DIP computed
    int dip_type;
    
  public:

    // @brief Learn one or two clauses in DIP learning
    bool learn_two_DIP_clauses;

  protected:
    ///////////////////////////////////////////////////////////////////////////////////////////
    // TEMPORARY VARIABLES
    //
    // these variables are allocated here to avoid repeated allocation costs

    /// @brief A learnt clause
    vec<Lit> mainLearnedClause;

    /// @brief Work stack for @code{litRedundant}: holds list of reason clauses to be examined
    vec<CRef> workStack;

    /// @brief Marks whether a variable has already been seen
    vec<bool> seen;

    /// @brief Mostly for efficient LBD computation. 'seen2[i]' will indicate if decision level or variable 'i' has been seen.
    vec<uint64_t> seen2;

    // @brief Marks whether a variable has already been seen (DIP learning)
    vec<bool> seen3;
      
    /// @brief Simple counter for marking purpose with 'seen2'.
    uint64_t counter;

    /// @brief Stores variables whose values have been set in @code{seen} and need to be
    /// cleared. Currently only used by @code{litRedundant}
    vec<Var> toClear;

    // CollectFirstUIP data structures

    vec<int> pathCs;

    vec<int> var_iLevel_tmp;

    vec<Var> involved_vars;

    ERManager* erManager;
      
  public:
    ///////////////////////////////////////////////////////////////////////////////////////////
    // STATISTICS

    /// @brief The total number of literals in from first UIP learnt clauses
    uint64_t max_literals;

    /// @brief The total number of literals in learnt clauses after clause simplification and
    /// minimization
    uint64_t tot_literals;

    /// @brief Total time spent in DIP computation
    double time_DIP;

    // @bried Number of conflicts where some DIP has been found (not necessarily learned)
    int conflicts_with_dip;

    // @bried Number of conflicts where a dangerous DIP situation is detected and hence no DIP-learning is done
    int conflicts_with_dangerous_dip;

  public:
    ///////////////////////////////////////////////////////////////////////////////////////////
    // CONSTRUCTORS

    /**
     * @brief Construct a new ConflictAnalyzer object
     * 
     * @param s Reference to main solver object
     */
    ConflictAnalyzer(Solver& s);

    /**
     * @brief Destroy the ConflictAnalyzer object
     * 
     */
    ~ConflictAnalyzer() = default;

  public:
    ///////////////////////////////////////////////////////////////////////////////////////////
    // STATE MODIFICATION

    /**
     * @brief Set up internal data structures for a new variable
     * 
     * @param v the variable to register
     */
    void newVar(Var v);

  public:
    ///////////////////////////////////////////////////////////////////////////////////////////
    // PUBLIC API

    /**
     * @brief Analyze the current conflict graph and generate a learnt clause and its
     * associated backtrack level. This uses DIP learning. If DIP learning is successful
     * out_learnt will contain the clause DIP -> conflict. This is the *main* learning clause,
     * as the BT point and the literal to be propagated is determined by this clause.
     * Additionally, a second clause out_learn_2 is generated that expresses UIP -> DIP.
     * Note that if no DIP is found (or we consider it is not good enough) the function 
     * returns false and in this case learnt_clause contains the 1UIP clause and out_learnt_UIP_to_DIP 
     * should be ignored

     * QUESTION: is out_learnt_UIP_to_DIP ordered in some sense?????

     * @param confl The clause responsible for the current conflict
     * @param out_learnt Output: the main learnt clause (DIP -> conflict)
     * @param out_btlevel Output: the backtrack level for the learnt clause
     * @param out_lbd Output: the LBD of the learnt clause
     * @param out_learnt2 Output: a secondary learnt clause (UIP -> DIP)
     * @return true if a proper DIP has been found
     * 
     * @pre 'out_learnt' is assumed to be cleared
     * @pre Current decision level must be greater than root level
     * 
     * @post 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
     * @post If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of
     * the rest of the literals. There may be others from the same level.
     */
    bool analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel, int& out_lbd, vec<Lit>& out_learnt_UIP_to_DIP);

    // Hack (Albert) We might want to make this nicer
    void notifyERManager(ERManager* erm);
      
    void simpleAnalyze(CRef confl, vec<Lit>& out_learnt, vec<CRef>& reason_clause, bool True_confl, int trailRecord);

    /**
     * @brief Specialized analysis procedure to express the final conflict in terms of
     * assumptions. Calculates the (possibly empty) set of assumptions that led to the
     * assignment of 'p'.
     * 
     * @param p the conflicting literal
     * @param out_conflict the set of assumptions leading to the assignment of @code{p}
     */
    void analyzeFinal(Lit p, vec<Lit>& out_conflict);
      
    bool collectFirstUIP(CRef confl);

  private:
    ///////////////////////////////////////////////////////////////////////////////////////////
    // HELPER FUNCTIONS

    /**
     * @brief Check whether a literal is redundant and can be removed.
     * 
     * @param p the literal to check for redundancy
     * @param abstract_levels an abstraction of decision levels, used to abort early if the
     * algorithm is visiting literals at levels that cannot be removed later.
     * @return true if p is redundant and can be removed, false otherwise
     * 
     * @note this is a helper method for @code{simplifyClauseDeep}
     * @pre the @code{workStack} vector is empty
     * @post the @code{workStack} vector is empty
     * 
     */
    bool litRedundant(Lit p, uint32_t abstract_levels);

    /**
     * @brief Simplify a learnt clause
     * 
     * @param simplified the output simplified clause.
     * @param toSimplify the clause to simplify.
     */
    void simplifyClauseDeep(vec<Lit>& simplified, const vec<Lit>& toSimplify);

    /**
     * @brief Check whether a reason clause is subsumed by a set of variables
     * 
     * @param c the reason clause for a variable in the learnt clause
     * @return true if the reason clause is subsumed, false otherwise
     * 
     * @note This is a helper method for @code{simplifyClauseBasic}
     * @note Precondition: @code{seen} is true for the variables in the learnt clause
     */
    bool reasonSubsumed(const Clause& c);

    /**
     * @brief Simplify a learnt clause
     * 
     * @param simplified the output simplified clause.
     * @param learntClause the learnt clause to modify.
     */
    void simplifyClauseBasic(vec<Lit>& simplified, const vec<Lit>& toSimplify);

    /**
     * @brief Learn a clause according to the first UIP learning scheme.
     * @param confl the conflict clause
     * @param learntClause the output learnt clause. Assumed to be empty initially.
     * @param learntClause2 the DIP->1UIP additional clause
     * @return true iff a DIP-learning has been performed
     */
    bool getDIPLearntClauses(CRef confl, vec<Lit>& learntClause, vec<Lit>& learntClause2);
    /**
     * @brief Computes the LBD of the clause DIP -> conflict
     * @param confl the conflict clause
     * @param x first member of the DIP
     * @param x second member of the DIP
     * @return LBD
     */
    int computeLBD_DIP2Conflict (CRef confl, Lit x, Lit y);
    
    /**
     * @brief Simplify a learnt clause
     * 
     * @param simplified the output simplified clause.
     * @param toSimplify the learnt clause to modify.
     */
    void simplifyClause(vec<Lit>& simplified, const vec<Lit>& toSimplify);

    /**
     * @brief Further learnt clause minimization by binary resolution.
     * 
     * @param out_learnt the learnt clause to modify.
     * @return true iff the learnt clause was modified.
     */
    bool binResMinimize(vec<Lit>& out_learnt);

    /**
     * @brief Enforce the watcher invariant for a learnt clause. Assumes that the variable with
     * the highest decision level is at index 0. Moves the variable with the next-highest
     * decision level to index 1.
     * 
     * @param learntClause the learnt clause to modify.
     */
    void enforceWatcherInvariant(vec<Lit>& learntClause);


    void writeDIPComputationInfo (TwoVertexBottlenecks& info, DIPGraphEncoder& encoder,
				  const vector<int>& predecessors, const vector<Lit>& predecessorLits,
				  const vector<int>& predIndex, const vector<Lit>& literalsInAnalysis, bool founDIP,
				  const vector<int>& pathA, const vector<int>& pathB);
    void writeEdgeInGraph (ofstream& out, const Lit& orig, const Lit& dest, bool colored);


    bool computeBestMiddleDIP (const TwoVertexBottlenecks& info, const DIPGraphEncoder& encoder, const vector<int>& pathA, const vector<int>& pathB, Lit& x, Lit& y);
    bool computeRandomDIP (int a, int b, TwoVertexBottlenecks& info, const DIPGraphEncoder& encoder, Lit& x, Lit& y);
    bool computeClosestDIPToConflict (int a, int b, TwoVertexBottlenecks& info, const DIPGraphEncoder& encoder, Lit& x, Lit& y);

    bool computeDIPClauses (int a, int b, CRef confl, TwoVertexBottlenecks& info, const DIPGraphEncoder& encoder, vec<Lit>& clause_to_learn1, vec<Lit>& clause_to_learn2, Lit UIP, vector<int>& pathA, vector<int>& pathB);
    bool ok_DIP (Lit dip1, Lit dip2, Lit UIP, CRef confl);
    bool checkSeen3();
  public:
    bool checkSeen();
  };
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // IMPLEMENTATION OF INLINE METHODS

  inline void ConflictAnalyzer::newVar(Var v) {
    seen.push(0);
    seen2.push(0);
    seen3.push(0);
    pathCs.push(0);
    var_iLevel_tmp.push(0);
  }
}

#endif
