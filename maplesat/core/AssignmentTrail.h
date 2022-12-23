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

#ifndef Minisat_AssignmentTrail_h
#define Minisat_AssignmentTrail_h

#include "core/SolverTypes.h"
#include "core/VariableDatabase.h"

namespace Minisat {
    class Solver;

    class AssignmentTrail {
    protected:
        struct VarData { CRef reason; int level; };
        static inline VarData mkVarData(CRef cr, int l) { VarData d = {cr, l}; return d; }

        vec<Lit>     trail;     // Assignment stack; stores all assigments made in the order they were made.
        vec<int>     trail_lim; // Separator indices for different decision levels in 'trail'.
        vec<VarData> vardata;   // Stores reason and level for each variable.

    public:
        //////////////////
        // CONSTRUCTORS //
        //////////////////

        AssignmentTrail(Solver* s);
        ~AssignmentTrail() = default;

        ////////////////
        // PUBLIC API //
        ////////////////

        void newVar(Var v);
        void relocAll(ClauseAllocator& to);

        void newDecisionLevel ();                              // Begins a new decision level.
        void uncheckedEnqueue (Lit p, CRef from = CRef_Undef); // Enqueue a literal. Assumes value of literal is undefined.
        bool enqueue          (Lit p, CRef from = CRef_Undef); // Test if fact 'p' contradicts current state, enqueue otherwise.
        void cancelUntil      (int level);                     // Backtrack until a certain level.

        int      decisionLevel()      const; // Gives the current decisionlevel.
        uint32_t abstractLevel(Var x) const; // Used to represent an abstraction of sets of decision levels.
        CRef     reason       (Var x) const;
        int      level        (Var x) const;

        bool locked (const Clause& c) const; // Returns TRUE if a clause is a reason for some implication in the current state.
        void handleEventClauseDeleted(const Clause& c);

        int nAssigns() const; // The current number of assigned literals.
        int nRootAssigns() const; // The current number of literals assigned at the root level

        Lit operator[] (int i) const;
        int indexOfDecisionLevel(int l) const;

        double progressEstimate () const; // DELETE THIS ?? IT'S NOT VERY USEFUL ...

    private:
        VariableDatabase& variableDatabase;
        ClauseAllocator& ca;
        Solver* solver;
    };

    inline void AssignmentTrail::newVar(Var v) {
        vardata.push(mkVarData(CRef_Undef, 0));
        trail  .capacity(v + 1);
    }

    inline void     AssignmentTrail::newDecisionLevel()                 { trail_lim.push(trail.size()); }
    inline bool     AssignmentTrail::enqueue         (Lit p, CRef from) { return variableDatabase.value(p) != l_Undef ? variableDatabase.value(p) != l_False : (uncheckedEnqueue(p, from), true); }
    
    inline int      AssignmentTrail::decisionLevel   ()      const      { return trail_lim.size(); }
    inline uint32_t AssignmentTrail::abstractLevel   (Var x) const      { return 1 << (level(x) & 31); }
    inline CRef     AssignmentTrail::reason          (Var x) const      { return vardata[x].reason; }
    inline int      AssignmentTrail::level           (Var x) const      { return vardata[x].level; }
    
    inline bool AssignmentTrail::locked (const Clause& c) const { return variableDatabase.value(c[0]) == l_True && reason(var(c[0])) != CRef_Undef && ca.lea(reason(var(c[0]))) == &c; }
    inline void AssignmentTrail::handleEventClauseDeleted(const Clause& c) { if (locked(c)) vardata[var(c[0])].reason = CRef_Undef; }

    inline int AssignmentTrail::nAssigns  () const { return trail.size(); }
    inline int AssignmentTrail::nRootAssigns() const { return trail_lim.size() == 0 ? trail.size() : trail_lim[0]; }

    inline Lit AssignmentTrail::operator[] (int i) const { return trail[i]; }
    inline int AssignmentTrail::indexOfDecisionLevel(int l) const { return l == 0 ? 0 : trail_lim[l - 1]; }
}

#endif