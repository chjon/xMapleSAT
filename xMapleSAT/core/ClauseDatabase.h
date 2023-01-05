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
        vec<CRef> learnts;

        /// @brief The maximum number of learnt clauses before clause deletion, which is triggered
        /// when the number of learnt clauses exceeds this value.
        double maxNumLearnts;

        /// @brief Timer for increasing the maximum size of the clause database. (initially 100
        /// by default)
        double learntSizeLimitGrowthTimer;

        /// @brief Number of conflicts remaining before the maximum size of the clause database
        /// should be increased
        int learntSizeLimitGrowthTimerCounter;

    private:
        ///////////////////////////////////////////////////////////////////////////////////////////////
        // TEMPORARY VARIABLES

        /// @brief Used to keep track of variables when computing LBD
        vec<uint64_t> lbd_seen;

    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // PARAMETERS

        /// @brief Indicates whether possibly inefficient linear scan for satisfied clauses should
        /// be performed in 'removeSatisfied'.
        bool remove_satisfied;

        /// @brief The fraction of wasted memory allowed before a garbage collection is triggered.
        double garbage_frac;

#if !LBD_BASED_CLAUSE_DELETION
        /// @brief A threshold activity for deleting clauses 
        double extra_lim;

        /// @brief Amount by which to decay clause activities
        double clause_decay;

        /// @brief Amount by which to bump the next clause
        double cla_inc;
#endif

        /// @brief The exponential growth factor for @code{learntSizeLimitGrowthTimer}.
        /// (default 1.5)
        double learntSizeTimerGrowthFactor;

#if !RAPID_DELETION
        /// @brief The initial limit for learnt clauses as a factor of the number of original
        /// clauses. (default 1 / 3)
        double learntSizeLimitFactorInitial;

        /// @brief The limit for learnt clauses is multiplied with this factor each restart.
        /// (default 1.1)
        double learntSizeLimitGrowthFactor;
#endif

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

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATE MODIFICATION

        /**
         * @brief Initialize the clause database size limit and associated timers.
         * 
         */
        void init(void);

        /**
         * @brief Set up internal data structures for a new variable
         * 
         * @param v the variable to register
         */
        void newVar(Var v);

        /**
         * @brief Add a clause to a database
         * 
         * @param ps the clause to add
         * @param db the database to which to add the clause
         * @param learnt true if @code{ps} is a learnt clause; false otherwise
         * @return the CRef of the newly added clause; CRef_Undef if the clause
         * is unit
         */
        template <typename V>
        CRef addClause(vec<Lit>& ps, V& db, bool learnt);

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
         * @brief Remove satisfied clauses from clause database.
         * 
         */
        void removeSatisfied(void);

        /**
         * @brief Reduce the set of learnt clauses
         * 
         */
        void checkReduceDB(void);

        /**
         * @brief Detach and free a clause
         * 
         * @param cr the clause to remove
         */
        void removeClause(CRef cr);

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

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // EVENT HANDLERS

        /**
         * @brief Update data structures for deletion heuristics when a clause appears in the
         * conflict graph.
         * 
         * @param c the clause that appears in the conflict graph
         */
        void handleEventClauseInConflictGraph(Clause& c);

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

    #if ! LBD_BASED_CLAUSE_DELETION
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
    #endif

        ///////////////////////////////////////////////////////////////////////////////////////////
        // HELPER FUNCTIONS

        /**
         * @brief Determine whether a clause should be deleted. This is a helper function for
         * @code{checkReduceDB}.
         * 
         * @details Remove half of the learnt clauses, minus the clauses locked by the current
         * assignment. Locked clauses are clauses that are reason to some assignment. Binary
         * clauses are never removed.
         * 
         * @param cdb The clause database object from which to remove clauses.
         * @param c the clause to check
         * @param index the index of the clause in the database
         * @return true iff the clause should be deleted
         * 
         * @note this function is declared static so that it can be passed as a function pointer.
         */
        static bool shouldRemoveReduceDB(
            const ClauseDatabase& cdb,
            const Clause& c,
            const int index
        );

        /**
         * @brief Determine whether a clause should be deleted. This is a helper function for
         * @code{removeSatisfied}.
         * 
         * @param cdb The clause database object from which to remove clauses.
         * @param c the clause to check
         * @param index the index of the clause in the database
         * @return true iff the clause should be deleted
         * 
         * @note this function is declared static so that it can be passed as a function pointer.
         */
        static bool shouldRemoveSatisfied(
            const ClauseDatabase& cdb,
            const Clause& c,
            const int index
        );

        /**
         * @brief Update clause database size limit data structures upon learning a clause.
         * 
         * @param cr the learnt clause
         * 
         */
        void handleEventLearntClause(CRef cr);

        /**
         * @brief Perform preprocessing of clause database to prepare for
         * clause deletion. This is a helper function for @code{checkReduceDB}.
         * 
         */
        void preprocessReduceDB(void);

        /**
         * @brief Remove clauses from a database according to a predicate. This is a helper
         * function for @code{checkReduceDB} and @code{removeSatisified}.
         * 
         * @tparam DeletionPredicate a function pointer: takes a clause and its index in the
         * database as arguments. Returns true iff the clause should be removed.
         * @param db the database from which to remove clauses
         * @param shouldRemove a predicate for deciding whether a clause should be deleted.
         */
        template <typename DeletionPredicate>
        void reduceDBWithPredicate(vec<CRef>& db, DeletionPredicate& shouldRemove);

        /**
         * @brief Reallocate clauses to defragment memory
         * 
         */
        void garbageCollect(void);

        /**
         * @brief Run garbage collection if memory is sufficiently fragmented
         * 
         * @param gf the fraction of wasted memory required to trigger garbage
         * collection.
         */
        void checkGarbage(double gf);

        /**
         * @brief Relocate all clauses
         * 
         * @param to the ClauseAllocator to relocate to
         */
        void relocAll(ClauseAllocator& to);
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE METHODS

    //////////////
    // ACCESSORS

    inline int ClauseDatabase::nClauses() const { return clauses.size(); }
    inline int ClauseDatabase::nLearnts() const { return learnts.size(); }

    ///////////////////////
    // STATE MODIFICATION

    inline void ClauseDatabase::newVar(Var v) {
        lbd_seen.push(0);
    }

    inline void ClauseDatabase::init(void) {
        // Initialize database size limit for clause deletion
    #if RAPID_DELETION
        maxNumLearnts = 2000;
    #else
        maxNumLearnts = nClauses() * learntSizeLimitFactorInitial;
    #endif

        // Initialize database growth timer
        learntSizeLimitGrowthTimerCounter = static_cast<int>(learntSizeLimitGrowthTimer);
    }

    template <typename V, typename T>
    static inline void vectorPush(V& vector, T& val) {
        vector.push(val);
    }

    template <>
    inline void vectorPush<std::vector<CRef>, CRef>(std::vector<CRef>& vector, CRef& val) {
        vector.push_back(val);
    }

    template <typename V>
    inline CRef ClauseDatabase::addClause(vec<Lit>& ps, V& db, bool learnt) {
        assert(ps.size() > 0);

        // Don't add unit clauses to the database -- they should be added to the trail instead
        if (ps.size() == 1) return CRef_Undef;
        
        // Allocate clause
        CRef cr = ca.alloc(ps, learnt);
        unitPropagator.attachClause(cr);
        vectorPush(db, cr);

        // Update stats
        if (learnt) learnts_literals += ps.size();
        else        clauses_literals += ps.size();

        return cr;
    }

    inline CRef ClauseDatabase::addLearntClause(vec<Lit>& ps) {
        const CRef cr = addClause(ps, learnts, true);
        handleEventLearntClause(cr);
        return cr;
    }

    inline CRef ClauseDatabase::addInputClause(vec<Lit>& ps) {
        assert(ps.size() > 1);
        return addClause(ps, clauses, false);
    }

    inline void ClauseDatabase::removeSatisfied(void) {
        // Remove satisfied clauses:
        reduceDBWithPredicate(learnts, shouldRemoveSatisfied);
        if (remove_satisfied)
            reduceDBWithPredicate(clauses, shouldRemoveSatisfied);

        // Perform garbage collection if needed
        checkGarbage(garbage_frac);
    }

    inline void ClauseDatabase::checkReduceDB() {
        if (nLearnts() - assignmentTrail.nAssigns() < maxNumLearnts) return;

        // Prepare for clause deletion
        preprocessReduceDB();

        // Perform clause deletion
        reduceDBWithPredicate(learnts, shouldRemoveReduceDB);

        // Perform garbage collection if needed
        checkGarbage(garbage_frac);

    #if RAPID_DELETION
        maxNumLearnts += 500;
    #endif
    }

    //////////////////////
    // UTILITY FUNCTIONS

    template<class V>
    inline int ClauseDatabase::lbd (const V& clause) {
        lbd_calls++;
        int lbd = 0;
        for (int i = 0; i < clause.size(); i++) {
            int l = assignmentTrail.level(var(clause[i]));
            if (lbd_seen[l] != lbd_calls) {
                lbd++;
                lbd_seen[l] = lbd_calls;
            }
        }
        return lbd;
    }

    ///////////////////
    // EVENT HANDLERS

    inline void ClauseDatabase::handleEventClauseInConflictGraph(Clause& c) {
    #if LBD_BASED_CLAUSE_DELETION
        if (c.learnt() && c.activity() > 2)
            c.activity() = lbd(c);
    #else
        if (c.learnt())
            claBumpActivity(c);
    #endif
    }

    ///////////
    // OUTPUT

    inline void ClauseDatabase::toDimacs(const char* file){ vec<Lit> as; toDimacs(file, as); }
    inline void ClauseDatabase::toDimacs(const char* file, Lit p){ vec<Lit> as; as.push(p); toDimacs(file, as); }
    inline void ClauseDatabase::toDimacs(const char* file, Lit p, Lit q){ vec<Lit> as; as.push(p); as.push(q); toDimacs(file, as); }
    inline void ClauseDatabase::toDimacs(const char* file, Lit p, Lit q, Lit r){ vec<Lit> as; as.push(p); as.push(q); as.push(r); toDimacs(file, as); }

    //////////////////////////////////////////
    // DELETION HEURISTIC STATE MODIFICATION

#if ! LBD_BASED_CLAUSE_DELETION
    inline void Solver::claDecayActivity() {
        cla_inc *= (1 / clause_decay);
    }

    inline void Solver::claBumpActivity (Clause& c) {
        const double RESCALE_THRESHOLD = 1e20;
        if ((c.activity() += cla_inc) <= RESCALE_THRESHOLD) continue;

        // Rescale:
        for (int i = 0; i < learnts.size(); i++)
            ca[learnts[i]].activity() /= RESCALE_THRESHOLD;
        cla_inc /= RESCALE_THRESHOLD;
    }
#endif

    /////////////////////
    // HELPER FUNCTIONS

    inline bool ClauseDatabase::shouldRemoveReduceDB(
        const ClauseDatabase& cdb,
        const Clause& c,
        const int index
    ) {
        // Don't delete binary or locked clauses. From the rest, delete clauses from the first half
        // and clauses with activity smaller than 'extra_lim'.
    #if LBD_BASED_CLAUSE_DELETION
        return
            c.activity() > 2 &&
            !cdb.assignmentTrail.locked(c) &&
            index < cdb.learnts.size() / 2;
    #else
        return
            c.size() > 2 &&
            !cdb.assignmentTrail.locked(c) &&
            (index < cdb.learnts.size() / 2 || c.activity() < cdb.extra_lim);
    #endif
    }

    inline bool ClauseDatabase::shouldRemoveSatisfied(
        const ClauseDatabase& cdb,
        const Clause& c,
        const int index
    ) {
        return cdb.assignmentTrail.satisfied(c);
    }

    template <typename DeletionPredicate>
    inline void ClauseDatabase::reduceDBWithPredicate(
        vec<CRef>& db,
        DeletionPredicate& shouldRemove
    ) {
        int i, j;
        for (i = j = 0; i < db.size(); i++){
            if (shouldRemove(*this, ca[db[i]], i))
                removeClause(db[i]);
            else
                db[j++] = db[i];
        }
        db.shrink(i - j);
    }

    inline void ClauseDatabase::relocAll(ClauseAllocator& to) {
        for (int i = 0; i < learnts.size(); i++) ca.reloc(learnts[i], to);
        for (int i = 0; i < clauses.size(); i++) ca.reloc(clauses[i], to);
    }

    inline void ClauseDatabase::checkGarbage(double gf) {
        if (ca.wasted() > ca.size() * gf)
            garbageCollect();
    }
}

#endif