#include <stdio.h>
#include "catch.hpp"
#include "Util.h"

namespace Minisat {

std::tr1::unordered_set<Lit> mkLitSet(std::initializer_list<int> elements) {
    std::tr1::unordered_set<Lit> s;
    for (auto element : elements) s.insert(mkLit(element));
    return s;
}

void setLitVec(vec<Lit>& v, std::initializer_list<int> elements) {
    v.clear();
    for (auto element : elements) v.push(mkLit(element));
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


}