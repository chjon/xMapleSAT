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
    SolverER::SolverER(Solver* s)
        : solver(s)
    {}

    SolverER::~SolverER() {}

    void SolverER::enforceWatcherInvariant(vec<Lit>& clause) {
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
            for (j = 1; j < clause.size(); j++)
                if (level(var(clause[j])) > level(var(clause[maxLvl_j])))
                    maxLvl_j = j;

            // Swap-in this literal at index 1:
            Lit tmp          = clause[1];
            clause[1]        = clause[maxLvl_j];
            clause[maxLvl_j] = tmp;
        }
    }

    void SolverER::introduceExtVars(std::tr1::unordered_map<Var, std::vector<CRef> >& ext_def_db) {
        if (m_extVarDefBuffer.size() == 0) return;

        extTimerStart();

        // Add extension variables
        // It is the responsibility of the user heuristic to ensure that we do not have pre-existing extension variables
        // for the provided literal pairs
        const bool ext_pref_sign = solver->ext_pref_sign;
        for (auto i = m_extVarDefBuffer.begin(); i != m_extVarDefBuffer.end(); i++) solver->newVar(ext_pref_sign);

        // Add extension definition clauses
        for (const ExtDef& def : m_extVarDefBuffer) {
            const Lit x = def.x, a = def.a, b = def.b;
            assert(var(x) > var(a) && var(x) > var(b));

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
        // for (ExtDef& def : m_extVarDefBuffer) er_prioritize(def.x);

        // Update stats and clean up
        total_ext_vars += m_extVarDefBuffer.size();
        max_ext_vars = std::max(max_ext_vars, total_ext_vars - deleted_ext_vars);
        m_extVarDefBuffer.clear();

        extTimerStop(ext_add_overhead);
    }

    void SolverER::addExtDefClause(std::vector<CRef>& db, Lit ext_lit, vec<Lit>& clause) {
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
            enforceWatcherInvariant(clause);

            // Add clause to data structures
            CRef cr = solver->ca.alloc(clause, false); // Allocating clause as if it were an original clause

            // Add clause to db
            db.push_back(cr);
            solver->attachClause(cr);
        }
    }
}