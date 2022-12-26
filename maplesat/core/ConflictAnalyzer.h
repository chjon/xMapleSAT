/****************************************************************************************[Solver.h]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010, Niklas Sorensson

Chanseok Oh's MiniSat Patch Series -- Copyright (c) 2015, Chanseok Oh
 
Maple_LCM, Based on MapleCOMSPS_DRUP -- Copyright (c) 2017, Mao Luo, Chu-Min LI, Fan Xiao: implementing a learnt clause minimisation approach
Reference: M. Luo, C.-M. Li, F. Xiao, F. Manya, and Z. L. , “An effective learnt clause minimization approach for cdcl sat solvers,” in IJCAI-2017, 2017, pp. to–appear.

Maple_LCM_Dist, Based on Maple_LCM -- Copyright (c) 2017, Fan Xiao, Chu-Min LI, Mao Luo: using a new branching heuristic called Distance at the beginning of search 

xMaple_LCM_Dist, based on Maple_LCM_Dist -- Copyright (c) 2022, Jonathan Chung, Vijay Ganesh, Sam Buss

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
#include "core/AssignmentTrail.h"
#include "core/BranchingHeuristicManager.h"

namespace Minisat {
    // Forward declarations
    class Solver;

    class ConflictAnalyzer {
    protected:
        ////////////////
        // PARAMETERS //
        ////////////////
        
        int ccmin_mode; // Controls conflict clause minimization (0=none, 1=basic, 2=deep).

        //////////////////////
        // MEMBER VARIABLES //
        //////////////////////

        vec<Lit> analyze_stack;

public:
        ////////////////
        // STATISTICS //
        ////////////////

        uint64_t max_literals;
        uint64_t tot_literals;

protected:
        ///////////////////////
        // SOLVER REFERENCES //
        ///////////////////////

        AssignmentTrail& assignmentTrail;
        BranchingHeuristicManager& branchingHeuristicManager;
        ClauseAllocator& ca;
        Solver* solver;

        //////////////////////
        // HELPER FUNCTIONS //
        //////////////////////

        bool litRedundant(Lit p, uint32_t abstract_levels); // (helper method for 'analyze()')

    public:
        //////////////////
        // CONSTRUCTORS //
        //////////////////

        ConflictAnalyzer(Solver* s);
        ~ConflictAnalyzer() = default;

        ////////////////
        // PUBLIC API //
        ////////////////

        void newVar(Var v);

        void analyze     (CRef confl, vec<Lit>& out_learnt, int& out_btlevel); // (bt = backtrack)
        void analyzeFinal(Lit p, vec<Lit>& out_conflict);                      // COULD THIS BE IMPLEMENTED BY THE ORDINARIY "analyze" BY SOME REASONABLE GENERALIZATION?
    };

    inline void ConflictAnalyzer::newVar(Var v) {
    }
}

#endif