#include "catch.hpp"

#include <stdio.h>
#include <initializer_list>
#include <core/SolverTypes.h>
#include <mtl/ExtDefMap.h>
#include <mtl/Vec.h>

namespace Minisat {

static std::tr1::unordered_set<Lit> mkLitSet(std::initializer_list<int> elements) {
    std::tr1::unordered_set<Lit> s;
    for (auto element : elements) s.insert(mkLit(element));
    return s;
}

static void setLitVec(vec<Lit>& v, std::initializer_list<int> elements) {
    v.clear();
    for (auto element : elements) v.push(mkLit(element));
}

TEST_CASE("Inserting extension variable definitions", "[ExtDefMap]") {
    Lit x = mkLit(0), a = mkLit(100), b = mkLit(200);
    ExtDefMap<Lit> xdm;
    REQUIRE(xdm.size() == 0);
    REQUIRE(xdm.degree(a) == 0);
    
    // No basis literals

    x = mkLit(1), a = mkLit(101), b = mkLit(201);
    xdm.insert(x, a, b);
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(x) == 0);
    REQUIRE(xdm.degree(a) == 1);
    REQUIRE(xdm.degree(b) == 1);

    x = mkLit(2), a = mkLit(102), b = mkLit(202);
    xdm.insert(x, a, b);
    REQUIRE(xdm.size() == 2);
    REQUIRE(xdm.degree(x) == 0);
    REQUIRE(xdm.degree(a) == 1);
    REQUIRE(xdm.degree(b) == 1);

    // One basis literal (any order)
    x = mkLit(3), a = mkLit(101), b = mkLit(203);
    xdm.insert(x, a, b);
    REQUIRE(xdm.degree(a) == 2);
    REQUIRE(xdm.degree(b) == 1);

    x = mkLit(4), a = mkLit(103), b = mkLit(201);
    xdm.insert(x, a, b);
    REQUIRE(xdm.degree(a) == 1);
    REQUIRE(xdm.degree(b) == 2);

    // Two basis literals (any order)
    x = mkLit(5), a = mkLit(101), b = mkLit(202);
    xdm.insert(x, a, b);
    REQUIRE(xdm.size() == 5);
    REQUIRE(xdm.degree(a) == 3);
    REQUIRE(xdm.degree(b) == 2);

    // Duplicate definitions should not change size or count
    x = mkLit(6), a = mkLit(101), b = mkLit(202);
    xdm.insert(x, a, b);
    REQUIRE(xdm.size() == 5);
    REQUIRE(xdm.degree(a) == 3);
    REQUIRE(xdm.degree(b) == 2);
}

TEST_CASE("Deleting extension variable definitions", "[ExtDefMap]") {
    Lit x = mkLit(0), a = mkLit(100), b = mkLit(200);
    Minisat::ExtDefMap<Lit> xdm;
    xdm.insert(x, a, b);

    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(a) == 1);
    REQUIRE(xdm.degree(b) == 1);

    // Zero basis literals
    xdm.erase(mkLit(101), mkLit(201));
    REQUIRE(xdm.size() == 1);

    // One basis literal (any order)
    xdm.erase(mkLit(100), mkLit(201));
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(a) == 1);
    REQUIRE(xdm.degree(b) == 1);

    xdm.erase(mkLit(101), mkLit(200));
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(a) == 1);
    REQUIRE(xdm.degree(b) == 1);

    ////////////////////////
    // Two basis literals //
    ////////////////////////

    xdm.insert(mkLit(1), mkLit(100), mkLit(201));
    xdm.insert(mkLit(2), mkLit(101), mkLit(200));
    xdm.insert(mkLit(3), mkLit(101), mkLit(202));
    REQUIRE(xdm.size() == 4);

    // Expect no change
    xdm.erase(mkLit(101), mkLit(201));
    REQUIRE(xdm.size() == 4);
    REQUIRE(xdm.degree(mkLit(101)) == 2);
    REQUIRE(xdm.degree(mkLit(201)) == 1);
    
    // Degrees should decrease by 1
    REQUIRE(xdm.degree(mkLit(100)) == 2);
    REQUIRE(xdm.degree(mkLit(200)) == 2);
    xdm.erase(mkLit(100), mkLit(200));
    REQUIRE(xdm.size() == 3);
    REQUIRE(xdm.degree(mkLit(100)) == 1);
    REQUIRE(xdm.degree(mkLit(200)) == 1);

    // No multiple deletion
    xdm.erase(mkLit(100), mkLit(200));
    REQUIRE(xdm.size() == 3);
    REQUIRE(xdm.degree(mkLit(100)) == 1);
    REQUIRE(xdm.degree(mkLit(200)) == 1);

    // Degrees should decrease by 1
    xdm.erase(mkLit(101), mkLit(200));
    REQUIRE(xdm.size() == 2);
    REQUIRE(xdm.degree(mkLit(101)) == 1);
    REQUIRE(xdm.degree(mkLit(200)) == 0);
    
    // Definition should be removed entirely
    xdm.erase(mkLit(101), mkLit(202));
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(mkLit(101)) == 0);
    REQUIRE(xdm.degree(mkLit(202)) == 0);

    // Should be able to delete in an order different from the definition
    xdm.erase(mkLit(201), mkLit(100));
    REQUIRE(xdm.size() == 0);
    REQUIRE(xdm.degree(mkLit(201)) == 0);
    REQUIRE(xdm.degree(mkLit(100)) == 0);
}

TEST_CASE("Clearing extension variable definitions", "[ExtDefMap]") {
    int x = 0;
    Minisat::ExtDefMap<Lit> xdm;
    xdm.insert(mkLit(x++), mkLit(100), mkLit(200));
    xdm.insert(mkLit(x++), mkLit(100), mkLit(201));
    xdm.insert(mkLit(x++), mkLit(101), mkLit(200));
    REQUIRE(xdm.size() == 3);

    xdm.clear();

    REQUIRE(xdm.size() == 0);
    REQUIRE(xdm.degree(mkLit(100)) == 0);
    REQUIRE(xdm.degree(mkLit(200)) == 0);
    REQUIRE(xdm.degree(mkLit(101)) == 0);
    REQUIRE(xdm.degree(mkLit(201)) == 0);
}

TEST_CASE("Deleting multiple extension variable definitions", "[ExtDefMap]") {
    Minisat::ExtDefMap<Lit> xdm;
    std::tr1::unordered_set<Lit> query;

    // Deleting a single extension variable with multiple definitions present
    xdm.insert(mkLit(0), mkLit(100), mkLit(200));
    xdm.insert(mkLit(1), mkLit(100), mkLit(201));
    xdm.insert(mkLit(2), mkLit(101), mkLit(200));
    xdm.erase(mkLitSet({0}));
    REQUIRE(xdm.size() == 2);
    REQUIRE(xdm.degree(mkLit(100)) == 1);
    REQUIRE(xdm.degree(mkLit(200)) == 1);
    REQUIRE(xdm.degree(mkLit(101)) == 1);
    REQUIRE(xdm.degree(mkLit(201)) == 1);

    // Reusing the same extension variable name
    xdm.insert(mkLit(0), mkLit(101), mkLit(201));
    xdm.erase(mkLitSet({0}));
    REQUIRE(xdm.size() == 2);
    REQUIRE(xdm.degree(mkLit(100)) == 1);
    REQUIRE(xdm.degree(mkLit(200)) == 1);
    REQUIRE(xdm.degree(mkLit(101)) == 1);
    REQUIRE(xdm.degree(mkLit(201)) == 1);

    // Deleting a single extension variable with no other definitions
    xdm.erase(mkLitSet({1}));
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(mkLit(100)) == 0);
    REQUIRE(xdm.degree(mkLit(200)) == 1);
    REQUIRE(xdm.degree(mkLit(101)) == 1);
    REQUIRE(xdm.degree(mkLit(201)) == 0);

    xdm.erase(mkLitSet({2}));
    REQUIRE(xdm.size() == 0);
    REQUIRE(xdm.degree(mkLit(200)) == 0);
    REQUIRE(xdm.degree(mkLit(101)) == 0);

    xdm.insert(mkLit(0), mkLit(100), mkLit(200));
    xdm.insert(mkLit(1), mkLit(100), mkLit(201));
    xdm.insert(mkLit(2), mkLit(101), mkLit(200));

    // Deleting the empty set
    xdm.erase(mkLitSet({}));
    REQUIRE(xdm.size() == 3);
    REQUIRE(xdm.degree(mkLit(100)) == 2);
    REQUIRE(xdm.degree(mkLit(200)) == 2);
    REQUIRE(xdm.degree(mkLit(101)) == 1);
    REQUIRE(xdm.degree(mkLit(201)) == 1);

    // Deleting multiple extension variables at once
    xdm.erase(mkLitSet({0, 1}));
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(mkLit(100)) == 0);
    REQUIRE(xdm.degree(mkLit(200)) == 1);
    REQUIRE(xdm.degree(mkLit(101)) == 1);
    REQUIRE(xdm.degree(mkLit(201)) == 0);

    // Deleting variables that are not in the map
    xdm.erase(mkLitSet({0, 1}));
    xdm.erase(mkLitSet({3}));
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(mkLit(100)) == 0);
    REQUIRE(xdm.degree(mkLit(200)) == 1);
    REQUIRE(xdm.degree(mkLit(101)) == 1);
    REQUIRE(xdm.degree(mkLit(201)) == 0);

    // Deleting where only some of the variables are in the map
    xdm.clear();
    xdm.insert(mkLit(0), mkLit(100), mkLit(200));
    xdm.insert(mkLit(1), mkLit(101), mkLit(201));
    xdm.insert(mkLit(2), mkLit(102), mkLit(202));
    xdm.erase(mkLitSet({0, 2, 3}));
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(mkLit(100)) == 0); REQUIRE(xdm.degree(mkLit(200)) == 0);
    REQUIRE(xdm.degree(mkLit(101)) == 1); REQUIRE(xdm.degree(mkLit(201)) == 1);
    REQUIRE(xdm.degree(mkLit(102)) == 0); REQUIRE(xdm.degree(mkLit(202)) == 0);

    // Deleting the last extension variable from the map
    xdm.erase(mkLitSet({0, 1, 2, 3}));
    REQUIRE(xdm.size() == 0);
    REQUIRE(xdm.degree(mkLit(100)) == 0); REQUIRE(xdm.degree(mkLit(200)) == 0);
    REQUIRE(xdm.degree(mkLit(101)) == 0); REQUIRE(xdm.degree(mkLit(201)) == 0);
    REQUIRE(xdm.degree(mkLit(102)) == 0); REQUIRE(xdm.degree(mkLit(202)) == 0);
}

TEST_CASE("Searching for extension variable definitions", "[ExtDefMap]") {
    Minisat::ExtDefMap<Lit> xdm;
    xdm.insert(mkLit(0), mkLit(100), mkLit(200));
    xdm.insert(mkLit(1), mkLit(100), mkLit(201));
    xdm.insert(mkLit(2), mkLit(101), mkLit(200));
    
    // Zero basis literals
    auto it = xdm.find(mkLit(102), mkLit(202));
    REQUIRE(it == xdm.end());

    // One basis literal (any order)
    it = xdm.find(mkLit(100), mkLit(202));
    REQUIRE(it == xdm.end());

    it = xdm.find(mkLit(102), mkLit(200));
    REQUIRE(it == xdm.end());

    // Two basis literals (any order)

    // No corresponding definition
    it = xdm.find(mkLit(101), mkLit(201));
    REQUIRE(it == xdm.end());

    // Existing corresponding definition
    it = xdm.find(mkLit(100), mkLit(201));
    REQUIRE(it != xdm.end());
    REQUIRE(it->second == mkLit(1));

    it = xdm.find(mkLit(200), mkLit(101));
    REQUIRE(it != xdm.end());
    REQUIRE(it->second == mkLit(2));
}

static void printVec(Minisat::vec<Lit>& v) {
    printf("[");
    if (v.size()) printf("%s%d", sign(v[0]) ? "-" : "", var(v[0]));
    for (int i = 1; i < v.size(); i++) printf(", %s%d", sign(v[i]) ? "-" : "", var(v[i]));
    printf("]");
}

static bool requireVecEqual(Minisat::vec<Lit>& actual, Minisat::vec<Lit>& expect) {
    REQUIRE(actual.size() == expect.size());
    if (actual.size() == expect.size()) {
        unsigned int numDifferences = 0;
        for (int i = 0; i < actual.size(); i++) {
            if (actual[i] != expect[i])
                numDifferences++;
        }
        REQUIRE(numDifferences == 0);
        if (numDifferences == 0) return true;
    }

    printf("Actual: "); printVec(actual); printf("\n");
    printf("Expect: "); printVec(expect); printf("\n");
    return false;
}

TEST_CASE("Substituting into clauses", "[ExtDefMap]") {
    Minisat::ExtDefMap<Lit> xdm;
    xdm.insert(mkLit(0), mkLit(100), mkLit(200));
    xdm.insert(mkLit(1), mkLit(100), mkLit(201));
    xdm.insert(mkLit(2), mkLit(101), mkLit(200));
    xdm.insert(mkLit(3), mkLit(102), mkLit(202));

    Minisat::vec<Lit> clause;
    Minisat::vec<Lit> expect;

    // Zero basis literals
    setLitVec(clause, {301, 302, 303, 304, 305});
    setLitVec(expect, {301, 302, 303, 304, 305});
    xdm.substitute(clause);
    REQUIRE(requireVecEqual(clause, expect));

    // Basis literals with no corresponding extension variables
    setLitVec(clause, {101, 300, 301, 302});
    setLitVec(expect, {101, 300, 301, 302});
    xdm.substitute(clause);
    REQUIRE(requireVecEqual(clause, expect));

    // Pair of basis literals with a corresponding extension variable
    setLitVec(clause, {100, 200, 300, 301, 302});
    setLitVec(expect, {  0,      300, 301, 302});
    xdm.substitute(clause);
    REQUIRE(requireVecEqual(clause, expect));

    // Multiple pairs of basis literals with corresponding extension variables
    setLitVec(clause, {100, 300, 202, 301, 302, 102, 200, 303});
    setLitVec(expect, {  0, 300,   3, 301, 302,           303});
    xdm.substitute(clause);
    REQUIRE(requireVecEqual(clause, expect));
}

}