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

    // Multiple undefined variables: undefined literals should be moved to the front
    setLitVec(clause, {100, 200, 101, 102, 103});
    setLitVec(prefix, {101, 102, 103});
    ser.enforceWatcherInvariant(clause);
    REQUIRE(requireVecPrefix(clause, prefix));

    setLitVec(clause, {101, 100, 102, 103, 200});
    setLitVec(prefix, {101, 102, 103});
    ser.enforceWatcherInvariant(clause);
    REQUIRE(requireVecPrefix(clause, prefix));

    // Single undefined variable: undefined literal should be moved to index 1
    setLitVec(clause, {100, 200, 101});
    ser.enforceWatcherInvariant(clause);
    REQUIRE(clause[0] != mkLit(101));
    REQUIRE(clause[1] == mkLit(101));

    // Zero undefined variables: literal with highest level should be in index 0; literal with second-highest level should be in index 1
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

}