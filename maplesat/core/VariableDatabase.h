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
    /**
     * @brief This class keeps track of variable assignments
     * 
     */
    struct VariableDatabase {
    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // MEMBER VARIABLES

        /// @brief The current variable assignments
        vec<lbool> assigns;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CONSTRUCTORS

        /**
         * @brief Construct a new VariableDatabase object
         * 
         */
        VariableDatabase() = default;

        /**
         * @brief Destroy the VariableDatabase object
         * 
         */
        ~VariableDatabase() = default;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // ACCESSORS

        /**
         * @brief Get the current number of variables
         * 
         * @return the current number of variables
         */
        int nVars (void) const;

        /**
         * @brief Get the truth assignment of a variable
         * 
         * @param x the variable whose truth assignment should be returned
         * @return l_True, l_False, or l_Undef, depending on the variable assignment
         */
        lbool value (Var x) const;

        /**
         * @brief Get the truth assignment of a literal
         * 
         * @param x the literal whose truth assignment should be returned
         * @return l_True, l_False, or l_Undef, depending on the underlying variable assignment
         */
        lbool value (Lit p) const;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATE MODIFICATION

        /**
         * @brief Set up internal data structures for a new variable
         * 
         * @return the ID of the new variable
         */
        Var newVar(void);

        /**
         * @brief Set the value of a variable
         * 
         * @param x the variable whose value should be set
         * @param val the truth assignment for the variable
         */
        void setVar(Var x, lbool val);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // UTILITY FUNCTIONS

        /**
         * @brief Check whether a clause is satisfied under the current variable assignment
         * 
         * @param c the clause to check
         * @return true iff a literal in the clause is satisfied
         */
        bool satisfied(const Clause& c) const;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE METHODS

    //////////////
    // ACCESSORS

    inline int VariableDatabase::nVars () const {
        return assigns.size();
    }

    inline lbool VariableDatabase::value (Var x) const {
        return assigns[x];
    }

    inline lbool VariableDatabase::value (Lit p) const {
        return assigns[var(p)] ^ sign(p);
    }

    ///////////////////////
    // STATE MODIFICATION

    inline Var VariableDatabase::newVar() {
        const Var v = nVars();
        assigns.push(l_Undef);
        return v;
    }

    inline void VariableDatabase::setVar(Var x, lbool val) {
        assigns[x] = val;
    }

    //////////////////////
    // UTILITY FUNCTIONS

    inline bool VariableDatabase::satisfied(const Clause& c) const {
        for (int i = 0; i < c.size(); i++)
            if (value(c[i]) == l_True) return true;
        return false;
    }
}

#endif