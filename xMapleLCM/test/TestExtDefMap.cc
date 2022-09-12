#include "catch.hpp"
#include <test/Util.h>
#include <core/SolverTypes.h>
#include <mtl/ExtDefMap.h>
#include <mtl/Vec.h>

namespace Minisat {

SCENARIO("Inserting extension variable definitions", "[ExtDefMap]") {
    GIVEN("An empty ExtDefMap") {
        ExtDefMap<Lit> xdm;
        REQUIRE(xdm.size() == 0);

        WHEN("adding a new definition") {
            Lit x = mkLit(1), a = mkLit(101), b = mkLit(201);
            xdm.insert(x, a, b);

            THEN("the size and degree increase") {
                REQUIRE(xdm.size() == 1);
                REQUIRE(xdm.degree(x) == 0);
                REQUIRE(xdm.degree(a) == 1);
                REQUIRE(xdm.degree(b) == 1);
            }
        }
    }

    GIVEN("An ExtDefMap with one definition") {
        ExtDefMap<Lit> xdm;
        Lit x1 = mkLit(1), a1 = mkLit(101), b1 = mkLit(201);
        xdm.insert(x1, a1, b1);
        REQUIRE(xdm.size() == 1);

        WHEN("adding a non-overlapping definition") {
            Lit x2 = mkLit(2), a2 = mkLit(102), b2 = mkLit(202);
            xdm.insert(x2, a2, b2);

            THEN("the size and degree increase only for the new definition") {
                REQUIRE(xdm.size() == 2);
                REQUIRE(xdm.degree(x1) == 0);
                REQUIRE(xdm.degree(a1) == 1);
                REQUIRE(xdm.degree(b1) == 1);
                REQUIRE(xdm.degree(x2) == 0);
                REQUIRE(xdm.degree(a2) == 1);
                REQUIRE(xdm.degree(b2) == 1);
            }
        }

        WHEN("adding a partially-overlapping definition") {
            Lit x2 = mkLit(2), a2 = mkLit(101), b2 = mkLit(202);
            xdm.insert(x2, a2, b2);

            THEN("the size and degree increase only for the new definition") {
                REQUIRE(xdm.size() == 2);
                REQUIRE(xdm.degree(x1) == 0);
                REQUIRE(xdm.degree(a1) == 2);
                REQUIRE(xdm.degree(b1) == 1);
                REQUIRE(xdm.degree(x2) == 0);
                REQUIRE(xdm.degree(b2) == 1);
            }
        }

        WHEN("adding a partially-overlapping definition (reversed order)") {
            Lit x2 = mkLit(2), a2 = mkLit(202), b2 = mkLit(101);
            xdm.insert(x2, a2, b2);

            THEN("the size and degree increase only for the new definition") {
                REQUIRE(xdm.size() == 2);
                REQUIRE(xdm.degree(x1) == 0);
                REQUIRE(xdm.degree(a1) == 2);
                REQUIRE(xdm.degree(b1) == 1);
                REQUIRE(xdm.degree(x2) == 0);
                REQUIRE(xdm.degree(a2) == 1);
            }
        }

        WHEN("adding a duplicate definition") {
            Lit x2 = mkLit(2), a2 = mkLit(101), b2 = mkLit(201);
            xdm.insert(x2, a2, b2);

            THEN("the size and degrees remain the same") {
                REQUIRE(xdm.size() == 1);
                REQUIRE(xdm.degree(x1) == 0);
                REQUIRE(xdm.degree(a1) == 1);
                REQUIRE(xdm.degree(b1) == 1);
                REQUIRE(xdm.degree(x2) == 0);
            }
        }

        WHEN("adding a duplicate definition (reversed order)") {
            Lit x2 = mkLit(2), a2 = mkLit(201), b2 = mkLit(101);
            xdm.insert(x2, a2, b2);

            THEN("the size and degrees remain the same") {
                REQUIRE(xdm.size() == 1);
                REQUIRE(xdm.degree(x1) == 0);
                REQUIRE(xdm.degree(a1) == 1);
                REQUIRE(xdm.degree(b1) == 1);
                REQUIRE(xdm.degree(x2) == 0);
            }
        }
    }

    GIVEN("An ExtDefMap with multiple definitions") {
        ExtDefMap<Lit> xdm;
        Lit x1 = mkLit(1), a1 = mkLit(101), b1 = mkLit(201); xdm.insert(x1, a1, b1);
        Lit x2 = mkLit(2), a2 = mkLit(102), b2 = mkLit(202); xdm.insert(x2, a2, b2);
        Lit x3 = mkLit(3), a3 = mkLit(101), b3 = mkLit(203); xdm.insert(x3, a3, b3);
        REQUIRE(xdm.size() == 3);

        WHEN("adding a non-overlapping definition") {
            Lit x = mkLit(4), a = mkLit(104), b = mkLit(204);
            xdm.insert(x, a, b);

            THEN("the size and degree increase only for the new definition") {
                REQUIRE(xdm.size() == 4);
                REQUIRE(xdm.degree(x) == 0);
                REQUIRE(xdm.degree(a) == 1);
                REQUIRE(xdm.degree(b) == 1);
                REQUIRE(xdm.degree(x1) == 0);
                REQUIRE(xdm.degree(x2) == 0);
                REQUIRE(xdm.degree(x3) == 0);
                REQUIRE(xdm.degree(a1) == 2);
                REQUIRE(xdm.degree(a2) == 1);
                REQUIRE(xdm.degree(b1) == 1);
                REQUIRE(xdm.degree(b2) == 1);
                REQUIRE(xdm.degree(b3) == 1);
            }
        }

        WHEN("adding a definition that overlaps with one other definition") {
            Lit x = mkLit(4), a = mkLit(102), b = mkLit(204);
            xdm.insert(x, a, b);

            THEN("the size and degree increase only for the new definition") {
                REQUIRE(xdm.size() == 4);
                REQUIRE(xdm.degree(x ) == 0);
                REQUIRE(xdm.degree(a ) == 2);
                REQUIRE(xdm.degree(b ) == 1);
                REQUIRE(xdm.degree(x1) == 0);
                REQUIRE(xdm.degree(x2) == 0);
                REQUIRE(xdm.degree(x3) == 0);
                REQUIRE(xdm.degree(a1) == 2);
                REQUIRE(xdm.degree(a2) == 2);
                REQUIRE(xdm.degree(b1) == 1);
                REQUIRE(xdm.degree(b2) == 1);
                REQUIRE(xdm.degree(b3) == 1);
            }
        }

        WHEN("adding a definition that overlaps with two other definitions") {
            Lit x = mkLit(4), a = mkLit(102), b = mkLit(203);
            xdm.insert(x, a, b);

            THEN("the size and degree increase only for the new definition") {
                REQUIRE(xdm.size() == 4);
                REQUIRE(xdm.degree(x ) == 0);
                REQUIRE(xdm.degree(a ) == 2);
                REQUIRE(xdm.degree(b ) == 2);
                REQUIRE(xdm.degree(x1) == 0);
                REQUIRE(xdm.degree(x2) == 0);
                REQUIRE(xdm.degree(x3) == 0);
                REQUIRE(xdm.degree(a1) == 2);
                REQUIRE(xdm.degree(a2) == 2);
                REQUIRE(xdm.degree(b1) == 1);
                REQUIRE(xdm.degree(b2) == 1);
            }
        }
    }
}

SCENARIO("Deleting extension variable definitions", "[ExtDefMap]") {
    GIVEN("an empty ExtDefMap") {
        Minisat::ExtDefMap<Lit> xdm;
        REQUIRE(xdm.size() == 0);

        WHEN("deleting a definition") {
            xdm.erase(mkLit(100), mkLit(200));

            THEN("the size and degrees do not change") {
                REQUIRE(xdm.size() == 0);
                REQUIRE(xdm.degree(mkLit(100)) == 0);
                REQUIRE(xdm.degree(mkLit(200)) == 0);
            }
        }
    }

    GIVEN("an ExtDefMap with one definition") {
        Minisat::ExtDefMap<Lit> xdm;
        Lit x1 = mkLit(0), a1 = mkLit(100), b1 = mkLit(200); xdm.insert(x1, a1, b1);
        REQUIRE(xdm.size() == 1);

        WHEN("deleting a non-existent definition with no overlap") {
            xdm.erase(mkLit(101), mkLit(201));

            THEN("the size and degrees do not change") {
                REQUIRE(xdm.size() == 1);
                REQUIRE(xdm.degree(mkLit(100)) == 1);
                REQUIRE(xdm.degree(mkLit(200)) == 1);
                REQUIRE(xdm.degree(mkLit(101)) == 0);
                REQUIRE(xdm.degree(mkLit(201)) == 0);
            }
        }

        WHEN("deleting a non-existent definition with partial (left) overlap") {
            xdm.erase(mkLit(100), mkLit(201));

            THEN("the size and degrees do not change") {
                REQUIRE(xdm.size() == 1);
                REQUIRE(xdm.degree(mkLit(100)) == 1);
                REQUIRE(xdm.degree(mkLit(200)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 0);
            }
        }

        WHEN("deleting a non-existent definition with partial (right) overlap") {
            xdm.erase(mkLit(101), mkLit(200));

            THEN("the size and degrees do not change") {
                REQUIRE(xdm.size() == 1);
                REQUIRE(xdm.degree(mkLit(100)) == 1);
                REQUIRE(xdm.degree(mkLit(200)) == 1);
                REQUIRE(xdm.degree(mkLit(101)) == 0);
            }
        }

        WHEN("deleting the definition") {
            xdm.erase(mkLit(100), mkLit(200));

            THEN("the size and degrees decrease by 1") {
                REQUIRE(xdm.size() == 0);
                REQUIRE(xdm.degree(mkLit(100)) == 0);
                REQUIRE(xdm.degree(mkLit(200)) == 0);
            }
        }

        WHEN("deleting the definition (reversed order)") {
            xdm.erase(mkLit(200), mkLit(100));

            THEN("the size and degrees decrease by 1") {
                REQUIRE(xdm.size() == 0);
                REQUIRE(xdm.degree(mkLit(100)) == 0);
                REQUIRE(xdm.degree(mkLit(200)) == 0);
            }
        }
    }

    GIVEN("an ExtDefMap with multiple definitions") {
        ExtDefMap<Lit> xdm;
        xdm.insert(mkLit(1), mkLit(101), mkLit(201));
        xdm.insert(mkLit(2), mkLit(102), mkLit(202));
        xdm.insert(mkLit(3), mkLit(101), mkLit(202));
        xdm.insert(mkLit(4), mkLit(104), mkLit(204));
        REQUIRE(xdm.size() == 4);

        WHEN("deleting a non-existent definition with no overlap") {
            xdm.erase(mkLit(105), mkLit(205));

            THEN("the size and degrees do not change") {
                REQUIRE(xdm.size() == 4);
                REQUIRE(xdm.degree(mkLit(101)) == 2);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 1);
                REQUIRE(xdm.degree(mkLit(202)) == 2);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }

        WHEN("deleting a non-existent definition with single overlap") {
            xdm.erase(mkLit(101), mkLit(205));

            THEN("the size and degrees do not change") {
                REQUIRE(xdm.size() == 4);
                REQUIRE(xdm.degree(mkLit(101)) == 2);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 1);
                REQUIRE(xdm.degree(mkLit(202)) == 2);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }

        WHEN("deleting a non-existent definition with double overlap") {
            xdm.erase(mkLit(101), mkLit(102));

            THEN("the size and degrees do not change") {
                REQUIRE(xdm.size() == 4);
                REQUIRE(xdm.degree(mkLit(101)) == 2);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 1);
                REQUIRE(xdm.degree(mkLit(202)) == 2);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }

        WHEN("deleting a definition with no overlap") {
            xdm.erase(mkLit(104), mkLit(204));

            THEN("the size and degrees decrease by 1") {
                REQUIRE(xdm.size() == 3);
                REQUIRE(xdm.degree(mkLit(101)) == 2);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 0);
                REQUIRE(xdm.degree(mkLit(201)) == 1);
                REQUIRE(xdm.degree(mkLit(202)) == 2);
                REQUIRE(xdm.degree(mkLit(204)) == 0);
            }
        }

        WHEN("deleting a definition with single (left) overlap") {
            xdm.erase(mkLit(101), mkLit(201));

            THEN("the size and degrees of the new definition decrease by 1") {
                REQUIRE(xdm.size() == 3);
                REQUIRE(xdm.degree(mkLit(101)) == 1);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 0);
                REQUIRE(xdm.degree(mkLit(202)) == 2);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }

        WHEN("deleting a definition with single (right) overlap") {
            xdm.erase(mkLit(102), mkLit(202));

            THEN("the size and degrees of the new definition decrease by 1") {
                REQUIRE(xdm.size() == 3);
                REQUIRE(xdm.degree(mkLit(101)) == 2);
                REQUIRE(xdm.degree(mkLit(102)) == 0);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 1);
                REQUIRE(xdm.degree(mkLit(202)) == 1);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }

        WHEN("deleting a definition with double overlap") {
            xdm.erase(mkLit(101), mkLit(202));

            THEN("the size and degrees of the new definition decrease by 1") {
                REQUIRE(xdm.size() == 3);
                REQUIRE(xdm.degree(mkLit(101)) == 1);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 1);
                REQUIRE(xdm.degree(mkLit(202)) == 1);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }

        WHEN("deleting in the reverse order of the definition") {
            xdm.erase(mkLit(101), mkLit(202));

            THEN("the size and degrees of the new definition decrease by 1") {
                REQUIRE(xdm.size() == 3);
                REQUIRE(xdm.degree(mkLit(101)) == 1);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 1);
                REQUIRE(xdm.degree(mkLit(202)) == 1);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }
    }
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

SCENARIO("Deleting multiple extension variable definitions", "[ExtDefMap]") {
    GIVEN("an ExtDefMap with multiple definitions") {
        Minisat::ExtDefMap<Lit> xdm;
        xdm.insert(mkLit(1), mkLit(101), mkLit(201));
        xdm.insert(mkLit(2), mkLit(102), mkLit(202));
        xdm.insert(mkLit(3), mkLit(101), mkLit(202));
        xdm.insert(mkLit(4), mkLit(104), mkLit(204));
        REQUIRE(xdm.size() == 4);

        WHEN("deleting the empty set") {
            xdm.erase(mkLitSet({}));

            THEN("the size and degrees do not change") {
                REQUIRE(xdm.size() == 4);
                REQUIRE(xdm.degree(mkLit(101)) == 2);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 1);
                REQUIRE(xdm.degree(mkLit(202)) == 2);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }

        WHEN("deleting a definition that does not exist") {
            xdm.erase(mkLitSet({5}));

            THEN("the size and degrees do not change") {
                REQUIRE(xdm.size() == 4);
                REQUIRE(xdm.degree(mkLit(101)) == 2);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 1);
                REQUIRE(xdm.degree(mkLit(202)) == 2);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }

        WHEN("deleting a single non-overlapping definition") {
            xdm.erase(mkLitSet({4}));

            THEN("the size and degrees decrease for the deleted definition") {
                REQUIRE(xdm.size() == 3);
                REQUIRE(xdm.degree(mkLit(101)) == 2);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 0);
                REQUIRE(xdm.degree(mkLit(201)) == 1);
                REQUIRE(xdm.degree(mkLit(202)) == 2);
                REQUIRE(xdm.degree(mkLit(204)) == 0);
            }
        }

        WHEN("deleting one single (left) overlapping definition") {
            xdm.erase(mkLitSet({1}));

            THEN("the size and degrees decrease for the deleted definition") {
                REQUIRE(xdm.size() == 3);
                REQUIRE(xdm.degree(mkLit(101)) == 1);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 0);
                REQUIRE(xdm.degree(mkLit(202)) == 2);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }

        WHEN("deleting one single (right) overlapping definition") {
            xdm.erase(mkLitSet({2}));

            THEN("the size and degrees decrease for the deleted definition") {
                REQUIRE(xdm.size() == 3);
                REQUIRE(xdm.degree(mkLit(101)) == 2);
                REQUIRE(xdm.degree(mkLit(102)) == 0);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 1);
                REQUIRE(xdm.degree(mkLit(202)) == 1);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }

        WHEN("deleting one double overlapping definition") {
            xdm.erase(mkLitSet({3}));

            THEN("the size and degrees decrease for the deleted definition") {
                REQUIRE(xdm.size() == 3);
                REQUIRE(xdm.degree(mkLit(101)) == 1);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 1);
                REQUIRE(xdm.degree(mkLit(202)) == 1);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }

        WHEN("deleting multiple non-overlapping definitions") {
            xdm.erase(mkLitSet({3, 4}));

            THEN("the size and degrees decrease for the deleted definitions") {
                REQUIRE(xdm.size() == 2);
                REQUIRE(xdm.degree(mkLit(101)) == 1);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 0);
                REQUIRE(xdm.degree(mkLit(201)) == 1);
                REQUIRE(xdm.degree(mkLit(202)) == 1);
                REQUIRE(xdm.degree(mkLit(204)) == 0);
            }
        }

        WHEN("deleting multiple overlapping definitions") {
            xdm.erase(mkLitSet({1, 3}));

            THEN("the size and degrees decrease for the deleted definitions") {
                REQUIRE(xdm.size() == 2);
                REQUIRE(xdm.degree(mkLit(101)) == 0);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 0);
                REQUIRE(xdm.degree(mkLit(202)) == 1);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }

        WHEN("deleting multiple overlapping definitions where only some variables are in the map") {
            xdm.erase(mkLitSet({1, 3, 5}));

            THEN("the size and degrees decrease for the deleted definitions") {
                REQUIRE(xdm.size() == 2);
                REQUIRE(xdm.degree(mkLit(101)) == 0);
                REQUIRE(xdm.degree(mkLit(102)) == 1);
                REQUIRE(xdm.degree(mkLit(104)) == 1);
                REQUIRE(xdm.degree(mkLit(201)) == 0);
                REQUIRE(xdm.degree(mkLit(202)) == 1);
                REQUIRE(xdm.degree(mkLit(204)) == 1);
            }
        }

        WHEN("deleting all the definitions from the map") {
            xdm.erase(mkLitSet({1, 2, 3, 4}));

            THEN("the size and degrees get set to 0") {
                REQUIRE(xdm.size() == 0);
                REQUIRE(xdm.degree(mkLit(101)) == 0);
                REQUIRE(xdm.degree(mkLit(102)) == 0);
                REQUIRE(xdm.degree(mkLit(104)) == 0);
                REQUIRE(xdm.degree(mkLit(201)) == 0);
                REQUIRE(xdm.degree(mkLit(202)) == 0);
                REQUIRE(xdm.degree(mkLit(204)) == 0);
            }
        }
    }
}

TEST_CASE("Searching for extension variable definitions", "[ExtDefMap]") {
    Minisat::ExtDefMap<Lit> xdm;
    xdm.insert(mkLit(1), mkLit(101), mkLit(201));
    xdm.insert(mkLit(2), mkLit(102), mkLit(202));
    xdm.insert(mkLit(3), mkLit(101), mkLit(202));
    xdm.insert(mkLit(4), mkLit(104), mkLit(204));
    REQUIRE(xdm.size() == 4);

    SECTION("searching for a non-existent definition with no overlap") {
        REQUIRE_FALSE(xdm.containsExt(mkLit(5)));
        REQUIRE_FALSE(xdm.containsPair(mkLit(105), mkLit(205)));
    }

    SECTION("searching for a non-existent definition with single overlap") {
        REQUIRE_FALSE(xdm.containsPair(mkLit(101), mkLit(200)));
        REQUIRE_FALSE(xdm.containsPair(mkLit(100), mkLit(201)));
    }

    SECTION("searching for a non-existent definition with double overlap") {
        REQUIRE_FALSE(xdm.containsPair(mkLit(102), mkLit(204)));
    }

    SECTION("searching for an existing definition") {
        REQUIRE(xdm.containsPair(mkLit(101), mkLit(201)));
        auto it = xdm.find(mkLit(101), mkLit(201));
        REQUIRE(it->second == mkLit(1));

        REQUIRE(xdm.containsPair(mkLit(201), mkLit(101)));
        it = xdm.find(mkLit(201), mkLit(101));
        REQUIRE(it->second == mkLit(1));
    }
}

SCENARIO("Substituting into clauses", "[ExtDefMap]") {
    GIVEN("an ExtDefMap with multiple definitions") {
        Minisat::vec<Lit> clause, extLits, expect1, expect2;
        Minisat::ExtDefMap<Lit> xdm;
        xdm.insert(mkLit(0), mkLit(100), mkLit(200));
        xdm.insert(mkLit(1), mkLit(100), mkLit(201));
        xdm.insert(mkLit(2), mkLit(101), mkLit(200));
        xdm.insert(mkLit(3), mkLit(102), mkLit(202));
    
        WHEN("substituting into a clause with no basis literals") {
            setLitVec(clause, {301, 302, 303, 304, 305});
            xdm.substitute(clause, extLits);

            THEN("the clause should be unchanged") {
                setLitVec(expect1, {301, 302, 303, 304, 305});
                setLitVec(expect2, {}); extLits.clear();
                REQUIRE(requireVecEqual(clause , expect1));
                REQUIRE(requireVecEqual(extLits, expect2));
            }
        }

        WHEN("substituting into a clause with basis literals but no corresponding extension definition") {
            setLitVec(clause, {101, 300, 301, 302});
            xdm.substitute(clause, extLits);

            THEN("the clause should be unchanged") {
                setLitVec(expect1, {101, 300, 301, 302});
                setLitVec(expect2, {});
                REQUIRE(requireVecEqual(clause , expect1));
                REQUIRE(requireVecEqual(extLits, expect2));
            }
        }

        WHEN("substituting into a clause containing a corresponding extension definition") {
            setLitVec(clause, {100, 200, 300, 301, 302});
            xdm.substitute(clause, extLits);

            THEN("the extension definition should be substituted into the clause") {
                setLitVec(expect1, {0, 300, 301, 302});
                setLitVec(expect2, {0});
                REQUIRE(requireVecEqual(clause , expect1));
                REQUIRE(requireVecEqual(extLits, expect2));
            }
        }

        WHEN("substituting into a clause containing an extension variable AND its corresponding extension definition") {
            setLitVec(clause, {0, 100, 200, 300, 301, 302});
            xdm.substitute(clause, extLits);

            THEN("the extension definition should not be substituted into the clause") {
                setLitVec(expect1, {0, 100, 200, 300, 301, 302});
                setLitVec(expect2, {0});
                REQUIRE(requireVecEqual(clause , expect1));
                REQUIRE(requireVecEqual(extLits, expect2));
            }
        }

        WHEN("substituting into a clause containing a corresponding extension definition (reverse order)") {
            setLitVec(clause, {200, 100, 300, 301, 302});
            xdm.substitute(clause, extLits);

            THEN("the extension definition should be substituted into the clause") {
                setLitVec(expect1, {0, 300, 301, 302});
                setLitVec(expect2, {0});
                REQUIRE(requireVecEqual(clause , expect1));
                REQUIRE(requireVecEqual(extLits, expect2));
            }
        }

        WHEN("substituting into a clause containing multiple corresponding extension definitions") {
            setLitVec(clause, {100, 300, 202, 301, 302, 102, 200, 303});
            xdm.substitute(clause, extLits);

            THEN("the extension definitions should be substituted into the clause") {
                setLitVec(expect1, {0, 300, 3, 301, 302, 303});
                setLitVec(expect2, {0, 3});
                REQUIRE(requireVecEqual(clause , expect1));
                REQUIRE(requireVecEqual(extLits, expect2));
            }
        }
    }
}

SCENARIO("Removing redundant literals from clauses", "[ExtDefMap]") {
    GIVEN("an ExtDefMap with multiple definitions") {
        Minisat::vec<Lit> clause, expect;
        Minisat::ExtDefMap<Lit> xdm;
        xdm.insert(mkLit(0), mkLit(100), mkLit(200));
        xdm.insert(mkLit(1), mkLit(100), mkLit(201));
        xdm.insert(mkLit(2), mkLit(101), mkLit(200));
        xdm.insert(mkLit(3), mkLit(102), mkLit(202));

        WHEN("simplifying a clause containing zero extension literals") {
            setLitVec(clause, {301, 302, 303, 305});
            xdm.absorb(clause);

            THEN("the clause should be unchanged") {
                setLitVec(expect, {301, 302, 303, 305});
                REQUIRE(requireVecEqual(clause, expect));
            }
        }
        
        WHEN("simplifying a clause containing zero basis literals") {
            setLitVec(clause, {301, 302, 303, 0, 305});
            xdm.absorb(clause);

            THEN("the clause should be unchanged") {
                setLitVec(expect, {301, 302, 303, 0, 305});
                REQUIRE(requireVecEqual(clause, expect));
            }
        }

        WHEN("simplifying a clause containing zero redundant literals") {
            setLitVec(clause, {101, 302, 303, 0, 305});
            xdm.absorb(clause);

            THEN("the clause should be unchanged") {
                setLitVec(expect, {101, 302, 303, 0, 305});
                REQUIRE(requireVecEqual(clause, expect));
            }
        }

        WHEN("simplifying a clause containing a single positive redundancy") {
            setLitVec(clause, {100, 302, 303, 0, 305});
            xdm.absorb(clause);

            THEN("the redundant basis literal should be removed") {
                setLitVec(expect, {302, 303, 0, 305});
                REQUIRE(requireVecEqual(clause, expect));
            }
        }

        WHEN("simplifying a clause containing two positive redundancies for the same extension variable") {
            setLitVec(clause, {200, 0, 302, 303, 100, 305});
            xdm.absorb(clause);

            THEN("the redundant literals should be removed") {
                setLitVec(expect, {0, 302, 303, 305});
                REQUIRE(requireVecEqual(clause, expect));
            }
        }

        WHEN("simplifying a clause containing positive redundancies for multiple extension variables") {
            setLitVec(clause, {200, 0, 302, 3, 100, 305, 102});
            xdm.absorb(clause);

            THEN("the redundant literals should be removed") {
                setLitVec(expect, {0, 302, 3, 305});
                REQUIRE(requireVecEqual(clause, expect));
            }
        }

        WHEN("simplifying a clause containing a single negative redundancy") {
            setLitVec(clause, {-100, 302, 303, -1, 305});
            xdm.absorb(clause);

            THEN("the redundant extension literal should be removed") {
                setLitVec(expect, {-100, 302, 303, 305});
                REQUIRE(requireVecEqual(clause, expect));
            }
        }

        WHEN("simplifying a clause containing two negative redundancies for the same extension variable") {
            setLitVec(clause, {-100, 302, 303, -1, -201, 305});
            xdm.absorb(clause);

            THEN("the redundant extension literal should be removed") {
                setLitVec(expect, {-100, 302, 303, -201, 305});
                REQUIRE(requireVecEqual(clause, expect));
            }
        }

        WHEN("simplifying a clause containing negative redundancies for multiple extension variables") {
            setLitVec(clause, {-100, -1, 302, -3, -201, 305, -102});
            xdm.absorb(clause);

            THEN("the redundant literals should be removed") {
                setLitVec(expect, {-100, 302, -201, 305, -102});
                REQUIRE(requireVecEqual(clause, expect));
            }
        }
    }
}

}