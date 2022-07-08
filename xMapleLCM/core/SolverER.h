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
#include "core/Solver.h"
#include "core/SolverTypes.h"
#include "mtl/Vec.h"

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

    int addExtDefClause(ClauseAllocator& ca, std::vector<CRef>& db, Lit ext_lit, vec<Lit>& clause);
    
    // Make sure the first two literals are in the right order for the watchers
    // Condition: clause must have length > 1 -- a unary clause should be propagated directly and not learnt
    void enforceWatcherInvariant(vec<Lit>& clause);

protected:
    // // Update stats
    // void updateExtFracStat(vec<Lit>& clause) {
    //     int numExtVarsInClause = getNumExtVars(clause);
    //     double extFrac = numExtVarsInClause / (double) clause.size();
    //     extfrac_total += extFrac;
    // }

    Solver* solver;

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

}

#endif