/*******************************************************************************************[SolverER.h]
Copyright (c) 2022, Jonathan Chung

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

#ifndef Minisat_SolverER_h
#define Minisat_SolverER_h

#include <map>
#include <vector>
#include <utility>
#include <mtl/Vec.h>
#include <mtl/ExtDefMap.h>
#include <core/Solver.h>
#include <core/SolverTypes.h>

namespace Minisat {

class SolverER {
public:
    SolverER(Solver* s);
    ~SolverER();

    inline int   level(Var x) const;
    inline lbool value(Var x) const;
    inline lbool value(Lit p) const;

#ifdef TESTING
    inline void set_value(Var x, lbool v, int l);
#endif

    using ProtoClause = vec<Lit>;
    using ExtDef = std::pair< Lit, std::vector<ProtoClause> >;

    // Clause Selection
    void selectClauses(std::vector<CRef>& selectedClauses);

    // Extension Variable Definition
    void defineExtVars(std::vector<ExtDef>& extVarDefs, const std::vector<CRef>& selectedClauses);

    // Extension Variable Introduction
    void introduceExtVars(const std::vector<ExtDef>& extVarDefs);

    /**
     * @brief Adds an extension definition clause to the appropriate clause database
     * 
     * @param ca Clause allocator to register the clause with the solver
     * @param db The clause database to which to add the clause
     * @param ext_lit The extension variable corresponding to the given clause
     * @param clause The vector of literals to add as a clause
     * 
     * @note Condition: current level must be zero
     */
    void addExtDefClause(ClauseAllocator& ca, std::vector<CRef>& db, Lit ext_lit, vec<Lit>& clause);

    /**
     * @brief Ensures the first two literals are in the right order for the watchers
     * 
     * @param clause a vector of multiple literals
     * 
     * @note @code{clause} must have length > 1. Unary clauses should be propagated directly and not learnt
     */
    void enforceWatcherInvariant(vec<Lit>& clause);

    // Extension Variable Substitution

    /**
     * @brief Substitute extension variables into a clause
     * 
     * @param clause The vector of literals in which to substitute
     */
    inline void substitute(vec<Lit>& clause) const;

protected:
    // // Update stats
    // void updateExtFracStat(vec<Lit>& clause) {
    //     int numExtVarsInClause = getNumExtVars(clause);
    //     double extFrac = numExtVarsInClause / (double) clause.size();
    //     extfrac_total += extFrac;
    // }

    Solver* solver;
    ExtDefMap<Lit> xdm;

    std::vector<CRef> m_selectedClauses;
    std::vector<ExtDef> m_extVarDefs;

#ifdef TESTING
    std::map< Var, std::pair<lbool, int> > test_value;
#endif

};

#ifdef TESTING
void  SolverER::set_value(Var x, lbool v, int l) { test_value.insert(std::make_pair(x, std::make_pair(v, l))); }
int   SolverER::level(Var x) const { auto it = test_value.find(x); return (it == test_value.end()) ? (-1     ) : (it->second.second); }
lbool SolverER::value(Var x) const { auto it = test_value.find(x); return (it == test_value.end()) ? (l_Undef) : (it->second.first ); }
lbool SolverER::value(Lit p) const { return value(var(p)) ^ sign(p); }
#else
int   SolverER::level(Var x) const { return solver->level(x); }
lbool SolverER::value(Var x) const { return solver->value(x); }
lbool SolverER::value(Lit p) const { return solver->value(p); }
#endif

void SolverER::substitute(vec<Lit>& clause) const { xdm.substitute(clause); }

}

#endif