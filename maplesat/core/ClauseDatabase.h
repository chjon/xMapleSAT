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
#include "core/UnitPropagator.h"

namespace Minisat {
    class Solver;
    class AssignmentTrail;
    class BranchingHeuristicManager;

    class ClauseDatabase {
    protected:
        //////////////////////
        // MEMBER VARIABLES //
        //////////////////////

        vec<CRef> clauses; // List of problem clauses.
        vec<CRef> learnts; // List of learnt clauses.

        /// @brief Indicates whether possibly inefficient linear scan for satisfied clauses
        // should be performed in 'simplify'.
        bool remove_satisfied;

        /// @brief The fraction of wasted memory allowed before a garbage collection is triggered.
        double garbage_frac;

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

        /**
         * @brief Get the current number of original (input) clauses.
         * 
         * @return The current number of original (input) clauses.
         */
        int nClauses(void) const;

        /**
         * @brief Get the current number of learnt clauses.
         * 
         * @return The current number of learnt clauses.
         */
        int nLearnts(void) const;

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

        /**
         * @brief Add a clause to a database
         * 
         * @param ps the clause to add
         * @param db the database to which to add the clause
         * @param learnt true if @code{ps} is a learnt clause; false otherwise
         * @return the CRef of the newly added clause; CRef_Undef if the clause is unit
         */
        CRef addClause(vec<Lit>& ps, vec<CRef>& db, bool learnt);

        /**
         * @brief Detach and free a clause
         * 
         * @param cr the clause to remove
         */
        void removeClause(CRef cr);

        /**
         * @brief Reallocate clauses to defragment memory
         * 
         */
        void garbageCollect(void);

        /**
         * @brief Run garbage collection if memory is sufficiently fragmented
         * 
         * @param gf the fraction of wasted memory required to trigger garbage collection
         */
        void checkGarbage(double gf);

        /**
         * @brief Remove satisfied clauses from a clause database
         * 
         * @param db the database from which to remove satisfied clauses
         */
        void removeSatisfied(vec<CRef>& db);

    public:
        //////////////////
        // CONSTRUCTORS //
        //////////////////

        /**
         * @brief Construct a new ClauseDatabase object
         * 
         * @param s Pointer to main solver object - must not be nullptr
         */
        ClauseDatabase(Solver* s);

        /**
         * @brief Destroy the ClauseDatabase object
         * 
         */
        ~ClauseDatabase() = default;

        ////////////////
        // PUBLIC API //
        ////////////////

        /**
         * @brief Add a clause to the learnt clause database
         * 
         * @param ps the list of literals to add as a clause
         * @return The CRef of the clause; CRef_Undef if ps is a unit clause
         */
        CRef addLearntClause(vec<Lit>& ps);

        /**
         * @brief Add a clause to the input clause database
         * 
         * @param ps the list of literals to add as a clause (at least two literals)
         * @return The CRef of the clause
         */
        CRef addInputClause(vec<Lit>& ps);

        /**
         * @brief Relocate all clauses
         * 
         * @param to the ClauseAllocator to relocate to
         */
        void relocAll(ClauseAllocator& to);

        /**
         * @brief Remove satisfied clauses from clause database.
         * 
         */
        void removeSatisfied(void);

        /**
         * @brief Reduce the set of learnt clauses
         * 
         */
        void reduceDB(void);

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

    ////////////////
    // STATISTICS //
    ////////////////

    inline int ClauseDatabase::nClauses() const { return clauses.size(); }
    inline int ClauseDatabase::nLearnts() const { return learnts.size(); }

    //////////////////////
    // HELPER FUNCTIONS //
    //////////////////////

    inline CRef ClauseDatabase::addClause(vec<Lit>& ps, vec<CRef>& db, bool learnt) {
        assert(ps.size() > 0);

        // Don't add unit clauses to the database -- they should be added to the trail instead
        if (ps.size() == 1) return CRef_Undef;
        
        // Allocate clause
        CRef cr = ca.alloc(ps, learnt);
        unitPropagator.attachClause(cr);
        db.push(cr);

        // Update stats
        if (learnt) learnts_literals += ps.size();
        else        clauses_literals += ps.size();

        return cr;
    }

    ////////////////
    // PUBLIC API //
    ////////////////

    inline CRef ClauseDatabase::addInputClause(vec<Lit>& ps) {
        assert(ps.size() > 1);
        return addClause(ps, clauses, false);
    }

    inline CRef ClauseDatabase::addLearntClause(vec<Lit>& ps) {
        return addClause(ps, learnts, true);
    }

    inline void ClauseDatabase::relocAll(ClauseAllocator& to) {
        for (int i = 0; i < learnts.size(); i++) ca.reloc(learnts[i], to);
        for (int i = 0; i < clauses.size(); i++) ca.reloc(clauses[i], to);
    }

    inline void ClauseDatabase::checkGarbage(double gf) {
        if (ca.wasted() > ca.size() * gf)
            garbageCollect();
    }

    inline void ClauseDatabase::removeSatisfied(void) {
        // Remove satisfied clauses:
        removeSatisfied(learnts);
        if (remove_satisfied)
            removeSatisfied(clauses);
        checkGarbage(garbage_frac);
    }

    ////////////
    // OUTPUT //
    ////////////

    inline void ClauseDatabase::toDimacs(const char* file){ vec<Lit> as; toDimacs(file, as); }
    inline void ClauseDatabase::toDimacs(const char* file, Lit p){ vec<Lit> as; as.push(p); toDimacs(file, as); }
    inline void ClauseDatabase::toDimacs(const char* file, Lit p, Lit q){ vec<Lit> as; as.push(p); as.push(q); toDimacs(file, as); }
    inline void ClauseDatabase::toDimacs(const char* file, Lit p, Lit q, Lit r){ vec<Lit> as; as.push(p); as.push(q); as.push(r); toDimacs(file, as); }
}

#endif