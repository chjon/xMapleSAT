#include <stdio.h>
#include "catch.hpp"
#include "Util.h"

namespace Minisat {

std::tr1::unordered_set<Lit> mkLitSet(const std::initializer_list<int>& elements) {
    std::tr1::unordered_set<Lit> s;
    for (const auto element : elements) s.insert(mkLit(element < 0 ? -element : element, element < 0));
    return s;
}

void setLitVec(vec<Lit>& v, const std::initializer_list<int>& elements) {
    v.clear();
    for (const auto element : elements) v.push(mkLit(element < 0 ? -element : element, element < 0));
}

static void printVec(Minisat::vec<Lit>& v) {
    printf("[");
    if (v.size()) printf("%s%d", sign(v[0]) ? "-" : "", var(v[0]));
    for (int i = 1; i < v.size(); i++) printf(", %s%d", sign(v[i]) ? "-" : "", var(v[i]));
    printf("]");
}

bool requireVecEqual(Minisat::vec<Lit>& actual, Minisat::vec<Lit>& expect) {
    unsigned int numDifferences = 0;
    if (actual.size() == expect.size()) {
        for (int i = 0; i < actual.size(); i++) {
            if (actual[i] != expect[i])
                numDifferences++;
        }
        if (numDifferences == 0) return true;
    }

    printf("Actual: "); printVec(actual); printf("\n");
    printf("Expect: "); printVec(expect); printf("\n");

    REQUIRE(actual.size() == expect.size());
    REQUIRE(numDifferences == 0);
    
    return false;
}

bool requireVecPrefix(Minisat::vec<Lit>& actual, Minisat::vec<Lit>& prefix) {
    unsigned int numDifferences = 0;
    if (actual.size() >= prefix.size()) {
        for (int i = 0; i < prefix.size(); i++) {
            if (actual[i] != prefix[i])
                numDifferences++;
        }
        if (numDifferences == 0) return true;
    }

    printf("Actual: "); printVec(actual); printf("\n");
    printf("Prefix: "); printVec(prefix); printf("\n");

    REQUIRE(actual.size() >= prefix.size());
    REQUIRE(numDifferences == 0);
    
    return false;
}

bool requireClauseEqual(const Clause& actual, const std::initializer_list<Lit>& elements) {
    // Ensure clause sizes are equal
    const unsigned int sz = static_cast<unsigned int>(actual.size());
    REQUIRE(sz == elements.size());
    if (sz != elements.size()) return false;

    // Ensure clauses contain the same elements
    unsigned int i = 0, numDiffs = 0;
    for (Lit l : elements) {
        if (actual[i++] != l) numDiffs++;
    }
    REQUIRE(numDiffs == 0);
    if (numDiffs != 0) return false;

    return true;
}

}