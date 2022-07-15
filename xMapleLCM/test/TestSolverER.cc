#include "catch.hpp"
#include <test/Util.h>
#include <core/SolverER.h>
#include <core/SolverTypes.h>
#include <mtl/Vec.h>

namespace Minisat {

TEST_CASE("Enforce watcher invariants", "[SolverER]") {
    SolverER ser(nullptr);
    ser.set_value(100, l_False, 0);
    ser.set_value(200, l_False, 1);
    ser.set_value(300, l_False, 2);
    ser.set_value(400, l_False, 2);

    vec<Lit> clause, prefix, expect;

    // Multiple unassigned variables: unassigned literals should be moved to the front
    setLitVec(clause, {100, 200, 101, 102, 103});
    setLitVec(prefix, {101, 102, 103});
    ser.enforceWatcherInvariant(clause);
    REQUIRE(requireVecPrefix(clause, prefix));

    setLitVec(clause, {101, 100, 102, 103, 200});
    setLitVec(prefix, {101, 102, 103});
    ser.enforceWatcherInvariant(clause);
    REQUIRE(requireVecPrefix(clause, prefix));

    // Single unassigned variable: unassigned literal should be moved to index 1 and highest-level literal should be in index 0
    setLitVec(clause, {100, 200, 101});
    setLitVec(prefix, {101, 200});
    ser.enforceWatcherInvariant(clause);
    REQUIRE(requireVecPrefix(clause, prefix));

    // Zero unassigned variables: literal with highest level should be in index 0; literal with second-highest level should be in index 1
    setLitVec(clause, {100, 200, 300});
    setLitVec(prefix, {200, 300});
    ser.enforceWatcherInvariant(clause);
    REQUIRE(requireVecPrefix(clause, prefix));

    setLitVec(clause, {100, 300, 200});
    setLitVec(prefix, {200, 300});
    ser.enforceWatcherInvariant(clause);
    REQUIRE(requireVecPrefix(clause, prefix));

    setLitVec(clause, {200, 400, 100, 300});
    setLitVec(prefix, {300, 400});
    ser.enforceWatcherInvariant(clause);
    REQUIRE(requireVecPrefix(clause, prefix));
}

TEST_CASE("Introducing extension variables", "[SolverER]") {
    Lit x, a, b;
    Solver s;
    SolverER& ser = *(s.ser);
    std::tr1::unordered_map<Var, std::vector<CRef> > db;
    std::vector< std::vector<Lit> > additional;

    // Set up variables for testing
    ser.originalNumVars = 10;
    for (int i = 0; i < ser.originalNumVars; i++) { s.newVar(); }

    // Test whether definition clauses are defined correctly
    x = mkLit(10), a = mkLit(0), b = mkLit(1);
    additional.clear(); additional.push_back(std::vector<Lit>({x, a, b, mkLit(2)}));
    ser.addToExtDefBuffer(ExtDef{ x, a, b, additional });
    ser.introduceExtVars(db);

    // Test whether clauses were added to the database
    std::tr1::unordered_map<Var, std::vector<CRef> >::iterator it = db.find(10);
    REQUIRE(it != db.end());

    if (it != db.end()) {
        // Test whether we get the expected number of clauses
        std::vector<CRef>& clauses = it->second;
        REQUIRE(clauses.size() == 4);

        if (clauses.size() == 4) {
            // Test whether we get the expected definition clauses
            REQUIRE(requireClauseEqual(s.getClause(clauses[0]), { ~x,  a,  b }));
            REQUIRE(requireClauseEqual(s.getClause(clauses[1]), {  x, ~a }));
            REQUIRE(requireClauseEqual(s.getClause(clauses[2]), {  x, ~b }));

            // Test whether we get the expected additional clause
            REQUIRE(requireClauseEqual(s.getClause(clauses[3]), { x, a, b, mkLit(2) }));
        }
    }

    // Test whether the extension variable was added to the solver
    REQUIRE(s.nVars() == ser.originalNumVars + 1);

    // Test whether the extension variable definition was stored in the extension definition map
    REQUIRE(ser.isCurrentExtVar(var(x)));

    // Check whether the buffer has been cleared
    REQUIRE(ser.extDefBufferSize() == 0);
}

TEST_CASE("Testing for valid definition pairs", "[SolverER]") {
    Solver s;
    SolverER& ser = *(s.ser);
    std::tr1::unordered_map<Var, std::vector<CRef> > db;
    std::tr1::unordered_set< std::pair<Lit, Lit> > generatedPairs;

    // Set up variables for testing
    ser.originalNumVars = 10;
    for (int i = 0; i < ser.originalNumVars; i++) { s.newVar(); }

    // Ensure literal pair consists of different variables
    Lit a = mkLit(1), b = mkLit(2);
    REQUIRE_FALSE(ser.isValidDefPair(mkLit(1, false), mkLit(1, true), generatedPairs));
    REQUIRE(ser.isValidDefPair(a, b, generatedPairs));
    REQUIRE(ser.isValidDefPair(b, a, generatedPairs));
    
    // Ensure literals in pair are not set at level 0
    ser.set_value(1, l_True, 0); ser.set_value(2, l_True, 0);
    REQUIRE_FALSE(ser.isValidDefPair(a, b, generatedPairs));
    ser.set_value(1, l_True, 0); ser.set_value(2, l_True, 1);
    REQUIRE_FALSE(ser.isValidDefPair(a, b, generatedPairs));
    ser.set_value(1, l_True, 1); ser.set_value(2, l_True, 0);
    REQUIRE_FALSE(ser.isValidDefPair(a, b, generatedPairs));
    ser.set_value(1, l_True, 1); ser.set_value(2, l_True, 1);
    REQUIRE      (ser.isValidDefPair(a, b, generatedPairs));

    // Ensure literal pair has not already been added
    generatedPairs.insert(mkLitPair(a, b));
    REQUIRE_FALSE(ser.isValidDefPair(a, b, generatedPairs));
    REQUIRE_FALSE(ser.isValidDefPair(b, a, generatedPairs));
    generatedPairs.clear();
    
    ser.addToExtDefBuffer(ExtDef{ mkLit(10), a, b, std::vector< std::vector<Lit> >() });
    ser.introduceExtVars(db);
    REQUIRE_FALSE(ser.isValidDefPair(b, a, generatedPairs));
}
}