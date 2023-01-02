/******************************************************************************[VariableDatabase.h]
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

#ifndef Minisat_VariableDatabase_h
#define Minisat_VariableDatabase_h

#include "core/SolverTypes.h"

namespace Minisat {
    struct VariableDatabase {
    protected:
        //////////////////////
        // MEMBER VARIABLES //
        //////////////////////

        vec<lbool> assigns; // The current variable assignments.

    public:
        //////////////////
        // CONSTRUCTORS //
        //////////////////

        VariableDatabase() = default;
        ~VariableDatabase() = default;

        ////////////////
        // PUBLIC API //
        ////////////////

        int   nVars ()      const     ; // The current number of variables.
        lbool value (Var x) const     ; // Gets the truth assignment of a variable
        lbool value (Lit p) const     ; // Gets the truth assignment of a literal
        Var   newVar()                ; // Introduce a new variable
        void  setVar(Var x, lbool val); // Set the value of a variable

        /**
         * @brief Check whether a clause is satisfied under the current variable assignment
         * 
         * @param c the clause to check
         * @return true iff a literal in the clause is satisfied
         */
        bool satisfied(const Clause& c) const;

    protected:
#ifdef TESTING
        inline void set_value(Var x, lbool v, int l);
#endif
    };

#ifdef TESTING
    inline void  VariableDatabase::set_value(Var x, lbool v, int l) {
        auto it = test_value.find(x);
        if (it == test_value.end()) test_value.insert(std::make_pair(x, std::make_pair(v, l)));
        else it->second = std::make_pair(v, l);
    }

    inline lbool VariableDatabase::value(Var x) const { auto it = test_value.find(x); return (it == test_value.end()) ? (l_Undef) : (it->second.first ); }
    inline lbool VariableDatabase::value(Lit p) const { return value(var(p)) ^ sign(p); }
#else
    inline int   VariableDatabase::nVars ()      const { return assigns.size(); }
    inline lbool VariableDatabase::value (Var x) const { return assigns[x]; }
    inline lbool VariableDatabase::value (Lit p) const { return assigns[var(p)] ^ sign(p); }
    inline Var   VariableDatabase::newVar() {
        const Var v = nVars();
        assigns.push(l_Undef);
        return v;
    }
    inline void VariableDatabase::setVar(Var x, lbool val) { assigns[x] = val; }
#endif

    inline bool VariableDatabase::satisfied(const Clause& c) const {
        for (int i = 0; i < c.size(); i++)
            if (value(c[i]) == l_True) return true;
        return false;
    }
}

#endif