#include "catch.hpp"
#include <test/Util.h>
#include <core/SolverER.h>
#include <core/SolverTypes.h>
#include <mtl/Vec.h>

namespace Minisat {

SCENARIO("Enforce watcher invariants", "[SolverER]") {
    GIVEN("Variable assignments") {
        vec<Lit> clause, prefix;
        SolverER ser(nullptr);
        ser.set_value(100, l_False, 0);
        ser.set_value(200, l_False, 1);
        ser.set_value(300, l_False, 2);
        ser.set_value(400, l_False, 2);

        WHEN("reordering a clause with multiple unassigned variables (packed together)") {
            setLitVec(clause, {100, 200, 101, 102, 103});
            ser.enforceWatcherInvariant(clause);

            THEN("unassigned literals should be moved to the front") {
                setLitVec(prefix, {101, 102, 103});
                REQUIRE_THAT(clause, vecPrefix(prefix));
            }
        }

        WHEN("reordering a clause with multiple unassigned variables (spread apart)") {
            setLitVec(clause, {101, 100, 102, 103, 200});
            ser.enforceWatcherInvariant(clause);

            THEN("unassigned literals should be moved to the front") {
                setLitVec(prefix, {101, 102, 103});
                REQUIRE_THAT(clause, vecPrefix(prefix));
            }
        }

        WHEN("reordering a clause with a single unassigned variable (highest-level literal before unassigned)") {
            setLitVec(clause, {100, 200, 101});
            ser.enforceWatcherInvariant(clause);

            THEN("the unassigned literal should be moved to index 1 and the highest-level literal should be in index 0") {
                setLitVec(prefix, {101, 200});
                REQUIRE_THAT(clause, vecPrefix(prefix));
            }
        }

        WHEN("reordering a clause with a single unassigned variable (highest-level literal after unassigned)") {
            setLitVec(clause, {100, 101, 200});
            ser.enforceWatcherInvariant(clause);

            THEN("the unassigned literal should be moved to index 1 and the highest-level literal should be in index 0") {
                setLitVec(prefix, {101, 200});
                REQUIRE_THAT(clause, vecPrefix(prefix));
            }
        }

        WHEN("reordering a clause with zero unassigned literals (highest-level before second-highest)") {
            setLitVec(clause, {100, 200, 300});
            ser.enforceWatcherInvariant(clause);

            THEN("the literal with highest level should be in index 0 and the literal with the second-highest level should be in index 1") {
                setLitVec(prefix, {200, 300});
                REQUIRE_THAT(clause, vecPrefix(prefix));
            }
        }

        WHEN("reordering a clause with zero unassigned literals (highest-level after second-highest)") {
            setLitVec(clause, {100, 300, 200});
            ser.enforceWatcherInvariant(clause);

            THEN("the literal with highest level should be in index 0 and the literal with the second-highest level should be in index 1") {
                setLitVec(prefix, {200, 300});
                REQUIRE_THAT(clause, vecPrefix(prefix));
            }
        }

        WHEN("reordering a clause with zero unassigned literals (multiple at highest-level)") {
            setLitVec(clause, {100, 400, 200, 300});
            ser.enforceWatcherInvariant(clause);

            THEN("the literals with the highest level should be in indices 0 and 1") {
                setLitVec(prefix, {300, 400});
                REQUIRE_THAT(clause, vecPrefix(prefix));
            }
        }
    }
}

SCENARIO("Introducing extension variables", "[SolverER]") {
    GIVEN("basis variables") {
        Solver s;
        SolverER& ser = *(s.ser);
        std::tr1::unordered_map<Var, std::vector<CRef> > db;
        std::vector< std::vector<Lit> > additional;
        vec<Lit> clause, expect;

        // Set up variables for testing
        ser.originalNumVars = 10;
        for (int i = 0; i < ser.originalNumVars; i++) { s.newVar(); }

        WHEN("introducing a new definition") {
            Lit x = mkLit(10), a = mkLit(0), b = mkLit(1);
            additional.push_back(std::vector<Lit>({x, a, b, mkLit(2)}));
            ser.addToExtDefBuffer(ExtDef{ x, a, b, additional });
            ser.introduceExtVars(db);

            THEN("clauses should be added to the database and the buffer should be cleared") {
                auto it = db.find(10);
                REQUIRE(it != db.end());
                if (it != db.end()) {
                    // Test whether we get the expected clauses
                    std::vector<CRef>& clauses = it->second;
                    REQUIRE(clauses.size() == 4);
                    clause2Vec(clause, s.ca[clauses[0]]); setVec(expect, { ~x,  a,  b           }); CHECK_THAT(clause, vecEqualUnordered(expect));
                    clause2Vec(clause, s.ca[clauses[1]]); setVec(expect, {  x, ~a               }); CHECK_THAT(clause, vecEqualUnordered(expect));
                    clause2Vec(clause, s.ca[clauses[2]]); setVec(expect, {  x,     ~b           }); CHECK_THAT(clause, vecEqualUnordered(expect));
                    clause2Vec(clause, s.ca[clauses[3]]); setVec(expect, {  x,  a,  b, mkLit(2) }); CHECK_THAT(clause, vecEqualUnordered(expect));
                }

                // Test whether the extension variable was added to the solver
                REQUIRE(s.nVars() == ser.originalNumVars + 1);

                // Test whether the extension variable definition was stored in the extension definition map
                REQUIRE(ser.isCurrentExtVar(var(x)));

                // Check whether the buffer was cleared
                REQUIRE(ser.extDefBufferSize() == 0);
            }
        }
    }
}

SCENARIO("Testing for valid definition pairs", "[SolverER]") {
    Solver s;
    SolverER& ser = *(s.ser);
    std::tr1::unordered_map<Var, std::vector<CRef> > db;

    // Set up variables for testing
    ser.originalNumVars = 10;
    for (int i = 0; i < ser.originalNumVars; i++) { s.newVar(); }

    GIVEN("no pre-existing pairs") {
        std::tr1::unordered_set< std::pair<Lit, Lit> > generatedPairs;

        WHEN("the input pair is two of the same variable") {
            Lit a = mkLit(1);

            THEN("the pair is not valid") {
                CHECK_FALSE(ser.isValidDefPair( a,  a, generatedPairs));
                CHECK_FALSE(ser.isValidDefPair( a, ~a, generatedPairs));
                CHECK_FALSE(ser.isValidDefPair(~a,  a, generatedPairs));
                CHECK_FALSE(ser.isValidDefPair(~a, ~a, generatedPairs));
            }
        }

        WHEN("the input pair consists of different variables") {
            Lit a = mkLit(1, false), b = mkLit(2, false);

            THEN("the pair is valid") {
                CHECK(ser.isValidDefPair( a,  b, generatedPairs));
                CHECK(ser.isValidDefPair(~a,  b, generatedPairs));
                CHECK(ser.isValidDefPair(~a, ~b, generatedPairs));
                CHECK(ser.isValidDefPair( a, ~b, generatedPairs));
            }
        }

        WHEN("checking input pairs with different assignment levels") {
            ser.set_value(1, l_True, 0); ser.set_value(2, l_False, 0);
            ser.set_value(3, l_True, 1); ser.set_value(4, l_False, 1);

            THEN("pairs containing a literal set at level 0 should be rejected") {
                CHECK_FALSE(ser.isValidDefPair(mkLit(1), mkLit(2), generatedPairs));
                CHECK_FALSE(ser.isValidDefPair(mkLit(1), mkLit(3), generatedPairs));
                CHECK_FALSE(ser.isValidDefPair(mkLit(1), mkLit(4), generatedPairs));
                CHECK_FALSE(ser.isValidDefPair(mkLit(3), mkLit(2), generatedPairs));
                CHECK_FALSE(ser.isValidDefPair(mkLit(4), mkLit(2), generatedPairs));
                CHECK      (ser.isValidDefPair(mkLit(3), mkLit(4), generatedPairs));
            }
        }
    }

    GIVEN("some pre-existing pairs") {
        std::tr1::unordered_set< std::pair<Lit, Lit> > generatedPairs;
        generatedPairs.insert(mkLitPair(mkLit(1), mkLit(2)));
        generatedPairs.insert(mkLitPair(mkLit(3), mkLit(4)));
        generatedPairs.insert(mkLitPair(mkLit(1), mkLit(4)));

        WHEN("the input pair has no overlap with existing pairs") {
            Lit a = mkLit(5), b = mkLit(6);
            
            THEN("the pair should be accepted") {
                CHECK(ser.isValidDefPair(a, b, generatedPairs));
                CHECK(ser.isValidDefPair(b, a, generatedPairs));
            }
        }

        WHEN("the input pair has single overlap with existing pairs") {
            Lit a = mkLit(1), b = mkLit(6);
            
            THEN("the pair should be accepted") {
                CHECK(ser.isValidDefPair(a, b, generatedPairs));
                CHECK(ser.isValidDefPair(b, a, generatedPairs));
            }
        }

        WHEN("the input pair has double overlap with existing pairs") {
            Lit a = mkLit(1), b = mkLit(3);
            
            THEN("the pair should be accepted") {
                CHECK(ser.isValidDefPair(a, b, generatedPairs));
                CHECK(ser.isValidDefPair(b, a, generatedPairs));
            }
        }

        WHEN("the input pair has already been added") {
            Lit a = mkLit(1), b = mkLit(2);
            
            THEN("the pair should be rejected") {
                CHECK_FALSE(ser.isValidDefPair(a, b, generatedPairs));
                CHECK_FALSE(ser.isValidDefPair(b, a, generatedPairs));
            }
        }
    }
}

SCENARIO("Choosing extension variables to delete", "[SolverER]") {

    GIVEN("Extension variables") {
        using DeletionPredicate = std::function<bool(Var)>;
        Solver s;
        SolverER& ser = *(s.ser);
        std::vector< std::vector<Lit> > additional;
        std::tr1::unordered_set<Lit> varsToDelete;
        vec<Lit> actual, expect;

        // Set up variables for testing
        ser.originalNumVars = 10;
        for (int i = 0; i < ser.originalNumVars; i++) { s.newVar(); }

        ser.addToExtDefBuffer(ExtDef{ mkLit(10), mkLit(1), mkLit(2), additional });
        ser.addToExtDefBuffer(ExtDef{ mkLit(11), mkLit(3), mkLit(4), additional });
        ser.addToExtDefBuffer(ExtDef{ mkLit(12), mkLit(5), mkLit(6), additional });
        ser.addToExtDefBuffer(ExtDef{ mkLit(13), mkLit(1), mkLit(4), additional });
        ser.introduceExtVars(ser.extDefs);

        WHEN("there are no restrictions") {
            DeletionPredicate deletionPredicate = [](Var v){ return true; };
            ser.getExtVarsToDelete(varsToDelete, deletionPredicate);

            THEN("all extension variables are chosen for deletion") {
                setVec(actual, varsToDelete);
                setLitVec(expect, {10, 11, 12, 13});
                REQUIRE_THAT(actual, vecEqualUnordered(expect));
            }
        }

        WHEN("some extension variables are basis literals") {
            DeletionPredicate deletionPredicate = [](Var v){ return true; };

            ser.addToExtDefBuffer(ExtDef{ mkLit(14), mkLit( 3), mkLit(10), additional });
            ser.addToExtDefBuffer(ExtDef{ mkLit(15), mkLit(11), mkLit(12), additional });
            ser.introduceExtVars(ser.extDefs);
            ser.getExtVarsToDelete(varsToDelete, deletionPredicate);

            THEN("only non-basis literal extension variables are chosen for deletion") {
                setVec(actual, varsToDelete);
                setLitVec(expect, {13, 14, 15});
                REQUIRE_THAT(actual, vecEqualUnordered(expect));
            }
        }

        WHEN("the deletion predicate rejects every variable") {
            DeletionPredicate deletionPredicate = [](Var v){ return false; };
            ser.getExtVarsToDelete(varsToDelete, deletionPredicate);

            THEN("no extension variables are chosen for deletion") {
                setVec(actual, varsToDelete);
                setLitVec(expect, {});
                REQUIRE_THAT(actual, vecEqualUnordered(expect));
            }
        }

        WHEN("the deletion predicate rejects some variables") {
            DeletionPredicate deletionPredicate = [](Var v){ return v % 2 == 0; };
            ser.getExtVarsToDelete(varsToDelete, deletionPredicate);

            THEN("none of the rejected variables are selected") {
                setVec(actual, varsToDelete);
                setLitVec(expect, {10, 12});
                REQUIRE_THAT(actual, vecEqualUnordered(expect));
            }
        }
    }
}

// SolverER::deleteExtVars
// SCENARIO("Deleting extension variables", "[SolverER]") {}
}