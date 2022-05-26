#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <stdio.h>
#include <mtl/ExtDefMap.h>
#include <mtl/Vec.h>

TEST_CASE("Inserting extension variable definitions", "[ExtDefMap]") {
    int x = 0, a = 100, b = 200;
    Minisat::ExtDefMap<int> xdm;
    REQUIRE(xdm.size() == 0);
    REQUIRE(xdm.degree(a) == 0);
    
    // No basis literals

    x++, a = 101, b = 201;
    xdm.insert(x, a, b);
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(x) == 0);
    REQUIRE(xdm.degree(a) == 1);
    REQUIRE(xdm.degree(b) == 1);

    x++, a = 102, b = 202;
    xdm.insert(x, a, b);
    REQUIRE(xdm.size() == 2);
    REQUIRE(xdm.degree(x) == 0);
    REQUIRE(xdm.degree(a) == 1);
    REQUIRE(xdm.degree(b) == 1);

    // One basis literal (any order)
    x++, a = 101, b = 203;
    xdm.insert(x, a, b);
    REQUIRE(xdm.degree(a) == 2);
    REQUIRE(xdm.degree(b) == 1);

    x++, a = 103, b = 201;
    xdm.insert(x, a, b);
    REQUIRE(xdm.degree(a) == 1);
    REQUIRE(xdm.degree(b) == 2);

    // Two basis literals (any order)
    x++, a = 101, b = 202;
    xdm.insert(x, a, b);
    REQUIRE(xdm.size() == 5);
    REQUIRE(xdm.degree(a) == 3);
    REQUIRE(xdm.degree(b) == 2);

    // Duplicate definitions should not change size or count
    x++, a = 101, b = 202;
    xdm.insert(x, a, b);
    REQUIRE(xdm.size() == 5);
    REQUIRE(xdm.degree(a) == 3);
    REQUIRE(xdm.degree(b) == 2);
}

TEST_CASE("Deleting extension variable definitions", "[ExtDefMap]") {
    int x = 0, a = 100, b = 200;
    Minisat::ExtDefMap<int> xdm;
    xdm.insert(x, a, b);

    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(a) == 1);
    REQUIRE(xdm.degree(b) == 1);

    // Zero basis literals
    xdm.erase(101, 201);
    REQUIRE(xdm.size() == 1);

    // One basis literal (any order)
    xdm.erase(a, b + 1);
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(a) == 1);
    REQUIRE(xdm.degree(b) == 1);

    xdm.erase(a + 1, b);
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(a) == 1);
    REQUIRE(xdm.degree(b) == 1);

    ////////////////////////
    // Two basis literals //
    ////////////////////////

    xdm.insert(x++, 100, 201);
    xdm.insert(x++, 101, 200);
    xdm.insert(x++, 101, 202);
    REQUIRE(xdm.size() == 4);

    // Expect no change
    xdm.erase(101, 201);
    REQUIRE(xdm.size() == 4);
    REQUIRE(xdm.degree(101) == 2);
    REQUIRE(xdm.degree(201) == 1);
    
    // Degrees should decrease by 1
    REQUIRE(xdm.degree(100) == 2);
    REQUIRE(xdm.degree(200) == 2);
    xdm.erase(100, 200);
    REQUIRE(xdm.size() == 3);
    REQUIRE(xdm.degree(100) == 1);
    REQUIRE(xdm.degree(200) == 1);

    // No multiple deletion
    xdm.erase(100, 200);
    REQUIRE(xdm.size() == 3);
    REQUIRE(xdm.degree(100) == 1);
    REQUIRE(xdm.degree(200) == 1);

    // Degrees should decrease by 1
    xdm.erase(101, 200);
    REQUIRE(xdm.size() == 2);
    REQUIRE(xdm.degree(101) == 1);
    REQUIRE(xdm.degree(200) == 0);
    
    // Definition should be removed entirely
    xdm.erase(101, 202);
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(101) == 0);
    REQUIRE(xdm.degree(202) == 0);

    // Should be able to delete in an order different from the definition
    xdm.erase(201, 100);
    REQUIRE(xdm.size() == 0);
    REQUIRE(xdm.degree(201) == 0);
    REQUIRE(xdm.degree(100) == 0);
}

TEST_CASE("Clearing extension variable definitions", "[ExtDefMap]") {
    int x = 0;
    Minisat::ExtDefMap<int> xdm;
    xdm.insert(x++, 100, 200);
    xdm.insert(x++, 100, 201);
    xdm.insert(x++, 101, 200);
    REQUIRE(xdm.size() == 3);

    xdm.clear();

    REQUIRE(xdm.size() == 0);
    REQUIRE(xdm.degree(100) == 0);
    REQUIRE(xdm.degree(200) == 0);
    REQUIRE(xdm.degree(101) == 0);
    REQUIRE(xdm.degree(201) == 0);
}

TEST_CASE("Deleting multiple extension variable definitions", "[ExtDefMap]") {
    Minisat::ExtDefMap<int> xdm;
    std::tr1::unordered_set<int> query;

    // Deleting a single extension variable with multiple definitions present
    xdm.insert(0, 100, 200);
    xdm.insert(1, 100, 201);
    xdm.insert(2, 101, 200);
    query.clear(); query.insert(0);
    xdm.erase(query);
    REQUIRE(xdm.size() == 2);
    REQUIRE(xdm.degree(100) == 1);
    REQUIRE(xdm.degree(200) == 1);
    REQUIRE(xdm.degree(101) == 1);
    REQUIRE(xdm.degree(201) == 1);

    // Reusing the same extension variable name
    xdm.insert(0, 101, 201);
    query.clear(); query.insert(0);
    xdm.erase(query);
    REQUIRE(xdm.size() == 2);
    REQUIRE(xdm.degree(100) == 1);
    REQUIRE(xdm.degree(200) == 1);
    REQUIRE(xdm.degree(101) == 1);
    REQUIRE(xdm.degree(201) == 1);

    // Deleting a single extension variable with no other definitions
    query.clear(); query.insert(1);
    xdm.erase(query);
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(100) == 0);
    REQUIRE(xdm.degree(200) == 1);
    REQUIRE(xdm.degree(101) == 1);
    REQUIRE(xdm.degree(201) == 0);

    query.clear(); query.insert(2);
    xdm.erase(query);
    REQUIRE(xdm.size() == 0);
    REQUIRE(xdm.degree(200) == 0);
    REQUIRE(xdm.degree(101) == 0);

    xdm.insert(0, 100, 200);
    xdm.insert(1, 100, 201);
    xdm.insert(2, 101, 200);

    // Deleting the empty set
    query.clear();
    xdm.erase(query);
    REQUIRE(xdm.size() == 3);
    REQUIRE(xdm.degree(100) == 2);
    REQUIRE(xdm.degree(200) == 2);
    REQUIRE(xdm.degree(101) == 1);
    REQUIRE(xdm.degree(201) == 1);

    // Deleting multiple extension variables at once
    query.clear(); query.insert(0); query.insert(1);
    xdm.erase(query);
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(100) == 0);
    REQUIRE(xdm.degree(200) == 1);
    REQUIRE(xdm.degree(101) == 1);
    REQUIRE(xdm.degree(201) == 0);

    // Deleting variables that are not in the map
    query.clear(); query.insert(0); query.insert(1);
    xdm.erase(query);
    query.clear(); query.insert(3);
    xdm.erase(query);
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(100) == 0);
    REQUIRE(xdm.degree(200) == 1);
    REQUIRE(xdm.degree(101) == 1);
    REQUIRE(xdm.degree(201) == 0);

    // Deleting where only some of the variables are in the map
    xdm.clear();
    xdm.insert(0, 100, 200);
    xdm.insert(1, 101, 201);
    xdm.insert(2, 102, 202);
    query.clear(); query.insert(0); query.insert(2); query.insert(3);
    xdm.erase(query);
    REQUIRE(xdm.size() == 1);
    REQUIRE(xdm.degree(100) == 0); REQUIRE(xdm.degree(200) == 0);
    REQUIRE(xdm.degree(101) == 1); REQUIRE(xdm.degree(201) == 1);
    REQUIRE(xdm.degree(102) == 0); REQUIRE(xdm.degree(202) == 0);

    // Deleting the last extension variable from the map
    query.clear(); query.insert(0); query.insert(1); query.insert(2); query.insert(3);
    xdm.erase(query);
    REQUIRE(xdm.size() == 0);
    REQUIRE(xdm.degree(100) == 0); REQUIRE(xdm.degree(200) == 0);
    REQUIRE(xdm.degree(101) == 0); REQUIRE(xdm.degree(201) == 0);
    REQUIRE(xdm.degree(102) == 0); REQUIRE(xdm.degree(202) == 0);
}

TEST_CASE("Searching for extension variable definitions", "[ExtDefMap]") {
    int x = 0;
    Minisat::ExtDefMap<int> xdm;
    xdm.insert(x++, 100, 200);
    xdm.insert(x++, 100, 201);
    xdm.insert(x++, 101, 200);
    
    // Zero basis literals
    auto it = xdm.find(102, 202);
    REQUIRE(it == xdm.end());

    // One basis literal (any order)
    it = xdm.find(100, 202);
    REQUIRE(it == xdm.end());

    it = xdm.find(102, 200);
    REQUIRE(it == xdm.end());

    // Two basis literals (any order)

    // No corresponding definition
    it = xdm.find(101, 201);
    REQUIRE(it == xdm.end());

    // Existing corresponding definition
    it = xdm.find(100, 201);
    REQUIRE(it != xdm.end());
    REQUIRE(it->second == 1);

    it = xdm.find(200, 101);
    REQUIRE(it != xdm.end());
    REQUIRE(it->second == 2);
}

static void printVec(Minisat::vec<int>& v) {
    printf("[");
    if (v.size()) printf("%d", v[0]);
    for (int i = 1; i < v.size(); i++) printf(", %d", v[i]);
    printf("]");
}

static bool requireVecEqual(Minisat::vec<int>& actual, Minisat::vec<int>& expect) {
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
    Minisat::ExtDefMap<int> xdm;
    xdm.insert(0, 100, 200);
    xdm.insert(1, 100, 201);
    xdm.insert(2, 101, 200);
    xdm.insert(3, 102, 202);

    Minisat::vec<int> clause;
    Minisat::vec<int> expect;

    // Zero basis literals
    clause.clear(); clause.push(301); clause.push(302); clause.push(303); clause.push(304); clause.push(305);
    expect.clear(); expect.push(301); expect.push(302); expect.push(303); expect.push(304); expect.push(305);
    xdm.substitute(clause);
    REQUIRE(requireVecEqual(clause, expect));

    // Basis literals with no corresponding extension variables
    clause.clear(); clause.push(100); clause.push(101); clause.push(300); clause.push(301); clause.push(302);
    expect.clear(); expect.push(100); expect.push(101); expect.push(300); expect.push(301); expect.push(302);
    xdm.substitute(clause);
    REQUIRE(requireVecEqual(clause, expect));

    // Pair of basis literals with a corresponding extension variable
    clause.clear(); clause.push(100); clause.push(200); clause.push(300); clause.push(301); clause.push(302);
    expect.clear(); expect.push(  0);                   expect.push(300); expect.push(301); expect.push(302);
    xdm.substitute(clause);
    REQUIRE(requireVecEqual(clause, expect));

    // Multiple pairs of basis literals with corresponding extension variables
    clause.clear(); clause.push(100); clause.push(300); clause.push(202); clause.push(301); clause.push(302); clause.push(102); clause.push(200); clause.push(303);
    expect.clear(); expect.push(  0); expect.push(300); expect.push(  3); expect.push(301); expect.push(302);                                     expect.push(303);
    xdm.substitute(clause);
    REQUIRE(requireVecEqual(clause, expect));
}