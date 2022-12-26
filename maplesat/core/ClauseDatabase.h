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

#ifndef Minisat_ClauseDatabase_h
#define Minisat_ClauseDatabase_h

#include <stdio.h>
#include "core/SolverTypes.h"
#include "core/VariableDatabase.h"

namespace Minisat {
    class Solver;
    class AssignmentTrail;
    class UnitPropagator;
    class BranchingHeuristicManager;

    class ClauseDatabase {
    protected:
        //////////////////////
        // MEMBER VARIABLES //
        //////////////////////

        vec<CRef> clauses; // List of problem clauses.
        vec<CRef> learnts; // List of learnt clauses.

        bool remove_satisfied; // Indicates whether possibly inefficient linear scan for satisfied clauses should be performed in 'simplify'.
        double garbage_frac; // The fraction of wasted memory allowed before a garbage collection is triggered.

        /////////////////////////
        // TEMPORARY VARIABLES //
        /////////////////////////

        vec<Lit> add_tmp; // Used to to reduce allocation overhead when adding new clauses

    public:
        ////////////////
        // STATISTICS //
        ////////////////

        uint64_t clauses_literals;
        uint64_t learnts_literals;

    protected:
        ///////////////////////
        // SOLVER REFERENCES //
        ///////////////////////
        
        VariableDatabase& variableDatabase;
        ClauseAllocator& ca;
        AssignmentTrail& assignmentTrail;
        UnitPropagator& unitPropagator;
        BranchingHeuristicManager& branchingHeuristicManager;
        Solver* solver;


        //////////////////////
        // HELPER FUNCTIONS //
        //////////////////////

        void removeClause(CRef cr); // Detach and free a clause.

    public:
        //////////////////
        // CONSTRUCTORS //
        //////////////////

        ClauseDatabase(Solver* s);
        ~ClauseDatabase() = default;

        ////////////////
        // PUBLIC API //
        ////////////////

        // Adding clauses

        bool addClause (const vec<Lit>& ps);  // Add a clause to the solver. 
        bool addEmptyClause();                // Add the empty clause, making the solver contradictory.
        bool addClause (Lit p);               // Add a unit clause to the solver. 
        bool addClause (Lit p, Lit q);        // Add a binary clause to the solver. 
        bool addClause (Lit p, Lit q, Lit r); // Add a ternary clause to the solver.

        /**
         * @brief Add a clause to the learnt clause database
         * 
         * @param ps the list of literals to add as a clause
         * @return The CRef of the generated clause; CRef_Undef if ps is a unit clause
         */
        CRef addLearntClause(vec<Lit>& ps);

        /**
         * @brief Add a clause to the solver without making superflous internal copy.
         * 
         * @param ps The list of literals to add as a clause. Will be modified.
         * @param learnt true if the clause is a learnt clause
         * @return true 
         * @return false 
         */
        bool addClause_(vec<Lit>& ps, bool learnt = false);

        void relocAll(ClauseAllocator& to);

        // Memory managment:
        //
        virtual void garbageCollect();
        void checkGarbage(double gf);
        void checkGarbage();

        // Statistics

        int nClauses() const; // The current number of original clauses.
        int nLearnts() const; // The current number of learnt clauses.

        void reduceDB       ();              // Reduce the set of learnt clauses.

        /**
         * @brief Remove satisfied clauses from clause database.
         * 
         */
        void removeSatisfied();

        void removeSatisfied(vec<CRef>& cs); // Shrink 'cs' to contain only non-satisfied clauses.

        ////////////
        // OUTPUT //
        ////////////

        void toDimacs(FILE* f, const vec<Lit>& assumps);            // Write CNF to file in DIMACS-format.
        void toDimacs(const char *file, const vec<Lit>& assumps);
        void toDimacs(FILE* f, Clause& c, vec<Var>& map, Var& max);

        // Convenience versions of 'toDimacs()':
        void toDimacs(const char* file);
        void toDimacs(const char* file, Lit p);
        void toDimacs(const char* file, Lit p, Lit q);
        void toDimacs(const char* file, Lit p, Lit q, Lit r);
    };

    //////////////////////////////////////
    // IMPLEMENTATION OF INLINE METHODS //
    //////////////////////////////////////

    inline bool ClauseDatabase::addClause     (const vec<Lit>& ps)  { ps.copyTo(add_tmp); return addClause_(add_tmp); }
    inline bool ClauseDatabase::addEmptyClause()                    { add_tmp.clear(); return addClause_(add_tmp); }
    inline bool ClauseDatabase::addClause     (Lit p)               { add_tmp.clear(); add_tmp.push(p); return addClause_(add_tmp); }
    inline bool ClauseDatabase::addClause     (Lit p, Lit q)        { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); return addClause_(add_tmp); }
    inline bool ClauseDatabase::addClause     (Lit p, Lit q, Lit r) { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); add_tmp.push(r); return addClause_(add_tmp); }

    inline void ClauseDatabase::relocAll(ClauseAllocator& to) {
        for (int i = 0; i < learnts.size(); i++) ca.reloc(learnts[i], to);
        for (int i = 0; i < clauses.size(); i++) ca.reloc(clauses[i], to);
    }

    inline void ClauseDatabase::checkGarbage(void){ return checkGarbage(garbage_frac); }
    inline void ClauseDatabase::checkGarbage(double gf){
        if (ca.wasted() > ca.size() * gf)
            garbageCollect();
    }

    inline void ClauseDatabase::removeSatisfied() {
        // Remove satisfied clauses:
        removeSatisfied(learnts);
        if (remove_satisfied)
            removeSatisfied(clauses);
        checkGarbage();
    }

    inline int ClauseDatabase::nClauses() const { return clauses.size(); }
    inline int ClauseDatabase::nLearnts() const { return learnts.size(); }

    inline void ClauseDatabase::toDimacs(const char* file){ vec<Lit> as; toDimacs(file, as); }
    inline void ClauseDatabase::toDimacs(const char* file, Lit p){ vec<Lit> as; as.push(p); toDimacs(file, as); }
    inline void ClauseDatabase::toDimacs(const char* file, Lit p, Lit q){ vec<Lit> as; as.push(p); as.push(q); toDimacs(file, as); }
    inline void ClauseDatabase::toDimacs(const char* file, Lit p, Lit q, Lit r){ vec<Lit> as; as.push(p); as.push(q); as.push(r); toDimacs(file, as); }
}

#endif