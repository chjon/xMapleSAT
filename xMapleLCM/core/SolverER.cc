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

namespace Minisat {
    SolverER::SolverER(Solver* s)
        : solver(s)
    {}

    SolverER::~SolverER() {}

    void SolverER::enforceWatcherInvariant(vec<Lit>& clause) {
        if (clause.size() == 1) return;

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

    // Add clause to extension definition database
    int SolverER::addExtDefClause(ClauseAllocator& ca, std::vector<CRef>& db, Lit ext_lit, vec<Lit>& clause) {
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

        if (clause.size() == 1) {
            solver->uncheckedEnqueue(ext_lit);
        } else {
            // Make sure the first two literals are in the right order for the watchers
            // shiftUnassigned(clause);

            // Add clause to data structures
            CRef cr = ca.alloc(clause, false); // Allocating clause as if it were an original clause

            // Add clause to db
            db.push_back(cr);
            solver->attachClause(cr);
        }

        return 0;
    }
}