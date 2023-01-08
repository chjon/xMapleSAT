/********************************************************************************[ClauseDatabase.h]
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

#ifndef Minisat_ClauseDatabase_h
#define Minisat_ClauseDatabase_h

// Don't change the actual numbers: clause.mark is a 2-bit value
// note: a value of 1 corresponds to a deleted clause
#define LOCAL 0
#define TIER2 2
#define CORE  3

#include <stdio.h>
#include <vector>

#include "core/AssignmentTrail.h"
#include "core/UnitPropagator.h"

namespace Minisat {
    // Forward declarations
    class Solver;

    /**
     * @brief This class manages clauses and clause deletion.
     * 
     */
    class ClauseDatabase {
    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // SOLVER REFERENCES
        
        ClauseAllocator& ca;
        AssignmentTrail& assignmentTrail;
        UnitPropagator& unitPropagator;
        Solver& solver;

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // MEMBER VARIABLES

        /// @brief List of input problem clauses.
        vec<CRef> clauses;

        /// @brief List of learnt clauses.
        vec<CRef> learnts_core;

        /// @brief List of learnt clauses.
        vec<CRef> learnts_tier2;

        /// @brief List of learnt clauses.
        vec<CRef> learnts_local;

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////////
        // TEMPORARY VARIABLES

    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // PARAMETERS

        /// @brief Indicates whether possibly inefficient linear scan for satisfied clauses should
        /// be performed in 'removeSatisfied'.
        bool remove_satisfied;

        /// @brief The fraction of wasted memory allowed before a garbage collection is triggered.
        double garbage_frac;

        /// @brief Conflict trigger threshold for the next tier2 clause database reduction
        uint64_t next_T2_reduce;

        /// @brief Conflict trigger threshold for the next local clause database reduction 
        uint64_t next_L_reduce;

        /// @brief Amount by which to bump the next clause
        double cla_inc;

        /// @brief Amount by which to decay clause activities
        double clause_decay;

        int core_lbd_cut;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATISTICS

        /// @brief The current total number of literals in original (input) clauses 
        uint64_t clauses_literals;

        /// @brief The current total number of literals in learnt clauses 
        uint64_t learnts_literals;
        
        /// @brief The total number of calls to LBD
        uint64_t lbd_calls;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CONSTRUCTORS

        /**
         * @brief Construct a new ClauseDatabase object
         * 
         * @param s Reference to main solver object
         */
        ClauseDatabase(Solver& s);

        /**
         * @brief Destroy the ClauseDatabase object
         * 
         */
        ~ClauseDatabase() = default;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // ACCESSORS

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

        template <int M>
        vec<CRef>& getDB();

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATE MODIFICATION

        /**
         * @brief Set up internal data structures for a new variable
         * 
         * @param v the variable to register
         */
        void newVar(Var v);

        /**
         * @brief Add a clause to the learnt clause database
         * 
         * @param ps the list of literals to add as a clause
         * @param lbd the lbd of the learnt clause
         * @param conflicts the current total number of conflicts seen by the solver
         * @return The CRef of the clause; CRef_Undef if ps is a unit clause
         */
        CRef addLearntClause(vec<Lit>& ps, int lbd, uint64_t conflicts);

        /**
         * @brief Add a clause to the input clause database
         * 
         * @param ps the list of literals to add as a clause (at least two literals)
         * @return The CRef of the clause
         */
        CRef addInputClause(vec<Lit>& ps);

        /**
         * @brief Remove satisfied clauses from clause database.
         * 
         */
        void removeSatisfied(void);

        /**
         * @brief Reduce the set of learnt clauses
         * 
         * @param conflicts the current total number of conflicts seen by the solver
         */
        void checkReduceDB(uint64_t conflicts);

        /**
         * @brief Detach and free a clause
         * 
         * @param cr the clause to remove
         */
        void removeClause(CRef cr);

        /**
         * @brief Relocate all clauses
         * 
         * @param to the ClauseAllocator to relocate to
         */
        void relocAll(ClauseAllocator& to);

        bool upgradeToCore(CRef cr, int lbd);

        /**
         * @brief Run garbage collection if memory is sufficiently fragmented
         * 
         */
        void checkGarbage(void);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // UTILITY FUNCTIONS

        /**
         * @brief Compute the LBD of a clause
         * 
         * @param clause the clause for which to compute LBD
         * @return the LBD of the clause
         */
        template<class V>
        int lbd (const V& clause);

        /**
         * @brief Get a list of all the learnt clauses that satisfy a predicate
         * 
         * @param output The output list
         * @param predicate the predicate function: takes a CRef as input and returns true iff the
         * clause satisfies the predicate. 
         */
        template<typename V, typename P>
        void selectLearntClauses(V& output, P predicate);

        void attachClause(CRef cr);

        void detachClause(CRef cr, bool strict);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // EVENT HANDLERS

        /**
         * @brief Update data structures for deletion heuristics when the solver encounters a
         * conflict.
         * 
         * @param conflicts 
         */
        void handleEventConflicted(uint64_t conflicts);

        /**
         * @brief Update data structures for deletion heuristics when a clause appears in the
         * conflict graph.
         * 
         * @param c the clause that appears in the conflict graph
         * @param conflicts the current total number of conflicts seen by the solver
         */
        void handleEventClauseInConflictGraph(CRef cr, uint64_t conflicts);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // OUTPUT

        ////////////////////////////////////////
        // Write CNF to file in DIMACS-format.

        void toDimacs(FILE* f, const vec<Lit>& assumps);
        void toDimacs(const char *file, const vec<Lit>& assumps);
        void toDimacs(FILE* f, Clause& c, vec<Var>& map, Var& max);

        //////////////////////////////////////////
        // Convenience versions of 'toDimacs()':
        
        void toDimacs(const char* file);
        void toDimacs(const char* file, Lit p);
        void toDimacs(const char* file, Lit p, Lit q);
        void toDimacs(const char* file, Lit p, Lit q, Lit r);

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // DELETION HEURISTIC STATE MODIFICATION

        /**
         * @brief Decay all clauses with the specified factor. 
         * 
         * @details Implemented by increasing the 'bump' value instead.
         * 
         */
        void claDecayActivity();
        
        /**
         * @brief Increase a clause with the current 'bump' value.
         * 
         * @param c the clause whose activity should be bumped
         */
        void claBumpActivity(Clause& c);

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // HELPER FUNCTIONS

        void removeSatisfied(vec<CRef>& cs);

        void safeRemoveSatisfied(vec<CRef>& cs, unsigned valid_mark);

        void reduceDB(void);
        void reduceDB_Tier2(void);

        /**
         * @brief Perform preprocessing of clause database to prepare for
         * clause deletion. This is a helper function for @code{checkReduceDB}.
         * 
         */
        void preprocessReduceDB(void);

        /**
         * @brief Reallocate clauses to defragment memory
         * 
         */
        void garbageCollect(void);
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE METHODS

    //////////////
    // ACCESSORS

    inline int ClauseDatabase::nClauses() const {
        return clauses.size();
    }

    inline int ClauseDatabase::nLearnts() const {
        return
            learnts_core.size() +
            learnts_tier2.size() +
            learnts_local.size();
    }

    template <> inline vec<CRef>& ClauseDatabase::getDB<CORE>() { return learnts_core; }
    template <> inline vec<CRef>& ClauseDatabase::getDB<TIER2>() { return learnts_tier2; }
    template <> inline vec<CRef>& ClauseDatabase::getDB<LOCAL>() { return learnts_local; }

    ///////////////////////
    // STATE MODIFICATION

    inline void ClauseDatabase::newVar(Var v) {
    }

    template <typename V, typename T>
    static inline void vectorPush(V& vector, T& val) {
        vector.push(val);
    }

    template <>
    inline void vectorPush<std::vector<CRef>, CRef>(std::vector<CRef>& vector, CRef& val) {
        vector.push_back(val);
    }

    inline CRef ClauseDatabase::addLearntClause(vec<Lit>& ps, int lbd, uint64_t conflicts) {
        if (ps.size() <= 1) return CRef_Undef;

        // Use LBD to select the appropriate learnt clause database
        CRef cr = ca.alloc(ps, true);
        ca[cr].set_lbd(lbd);
        if (lbd <= core_lbd_cut) {
            learnts_core.push(cr);
            ca[cr].mark(CORE);
        } else if (lbd <= 6) {
            learnts_tier2.push(cr);
            ca[cr].mark(TIER2);
            ca[cr].touched() = conflicts;
        } else {
            learnts_local.push(cr);
            claBumpActivity(ca[cr]);
        }

        attachClause(cr);
        claDecayActivity();
        return cr;
    }

    inline CRef ClauseDatabase::addInputClause(vec<Lit>& ps) {
        if (ps.size() <= 1) return CRef_Undef;
        CRef cr = ca.alloc(ps, false);
        clauses.push(cr);
        attachClause(cr);
        return cr;
    }

    inline void ClauseDatabase::checkReduceDB(uint64_t conflicts) {
        if (conflicts >= next_T2_reduce) {
            next_T2_reduce = conflicts + 10000;
            reduceDB_Tier2();
        }
        if (conflicts >= next_L_reduce) {
            next_L_reduce = conflicts + 15000;
            reduceDB();
        }
    }

    inline void ClauseDatabase::relocAll(ClauseAllocator& to) {
        // All learnt:
        for (int i = 0; i < learnts_core.size(); i++)  ca.reloc(learnts_core[i], to);
        for (int i = 0; i < learnts_tier2.size(); i++) ca.reloc(learnts_tier2[i], to);
        for (int i = 0; i < learnts_local.size(); i++) ca.reloc(learnts_local[i], to);

        // All original:
        int i, j;
        for (i = j = 0; i < clauses.size(); i++) {
            if (ca[clauses[i]].mark() == 1) continue;

            ca.reloc(clauses[i], to);
            clauses[j++] = clauses[i];
        }
        clauses.shrink(i - j);
    }

    inline bool ClauseDatabase::upgradeToCore(CRef cr, int lbd) {
        if (lbd > core_lbd_cut) return false;
        learnts_core.push(cr);
        ca[cr].mark(CORE);
        return true;
    }

    inline void ClauseDatabase::checkGarbage(void) {
        if (ca.wasted() > ca.size() * garbage_frac)
            garbageCollect();
    }

    //////////////////////
    // UTILITY FUNCTIONS

    template<typename V, typename P>
    void ClauseDatabase::selectLearntClauses(V& output, P predicate) {

    }

    ///////////////////
    // EVENT HANDLERS

    inline void ClauseDatabase::handleEventConflicted(uint64_t conflicts) {
        if (conflicts == 100000 && learnts_core.size() < 100)
            core_lbd_cut = 5;
    }

    inline void ClauseDatabase::handleEventClauseInConflictGraph(CRef cr, uint64_t conflicts) {
        Clause& c = ca[cr];

        // Skip input and core clauses
        if (!c.learnt() || c.mark() == CORE) return;

        // Update LBD
        int lbd = assignmentTrail.computeLBD(c);
        if (lbd < c.lbd()) {
            if (c.lbd() <= 30)
                c.removable(false); // Protect once from reduction.
            
            c.set_lbd(lbd);
            if (lbd <= core_lbd_cut) {
                learnts_core.push(cr);
                c.mark(CORE);
            } else if (lbd <= 6 && c.mark() == LOCAL) {
                // Bug: 'cr' may already be in 'learnts_tier2', e.g., if 'cr' was demoted from TIER2
                // to LOCAL previously and if that 'cr' is not cleaned from 'learnts_tier2' yet.
                learnts_tier2.push(cr);
                c.mark(TIER2);
            }
        }

        if (c.mark() == TIER2) {
            c.touched() = conflicts;
        } else if (c.mark() == LOCAL) {
            claBumpActivity(c);
        }
    }

    ///////////
    // OUTPUT

    inline void ClauseDatabase::toDimacs(const char* file){ vec<Lit> as; toDimacs(file, as); }
    inline void ClauseDatabase::toDimacs(const char* file, Lit p){ vec<Lit> as; as.push(p); toDimacs(file, as); }
    inline void ClauseDatabase::toDimacs(const char* file, Lit p, Lit q){ vec<Lit> as; as.push(p); as.push(q); toDimacs(file, as); }
    inline void ClauseDatabase::toDimacs(const char* file, Lit p, Lit q, Lit r){ vec<Lit> as; as.push(p); as.push(q); as.push(r); toDimacs(file, as); }

    //////////////////////////////////////////
    // DELETION HEURISTIC STATE MODIFICATION

    inline void ClauseDatabase::claDecayActivity() {
        cla_inc *= (1 / clause_decay);
    }

    inline void ClauseDatabase::claBumpActivity (Clause& c) {
        const double RESCALE_THRESHOLD = 1e20;
        if ((c.activity() += cla_inc) <= RESCALE_THRESHOLD) return;

        // Rescale:
        for (int i = 0; i < learnts_local.size(); i++)
            ca[learnts_local[i]].activity() /= RESCALE_THRESHOLD;
        cla_inc /= RESCALE_THRESHOLD;
    }
}

#endif