#include "catch.hpp"
#include <test/Util.h>
#include <core/SolverER.h>
#include <core/SolverTypes.h>
#include <mtl/Vec.h>

namespace Minisat {

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

static void setVariables(SolverER& ser, int i_undef, int i_max, int numVars) {
    ser.test_value.clear();
    for (int i = 0; i < numVars; i++)
        if (i != i_undef && i != i_max)
            ser.set_value(i, l_False, i);
    ser.set_value(i_max, l_False, numVars);
}

SCENARIO("Find asserting extension definition clause", "[SolverER]") {
    GIVEN("Variables") {
        Solver s;
        SolverER& ser = *(s.ser);
        std::tr1::unordered_set<Lit> varsToDelete;
        vec<Lit> ps, actual, expect;
        int i_undef = -1, i_max = -1;
        CRef cr;

        // Set up variables for testing
        ser.originalNumVars = 10;
        for (int i = 0; i < ser.originalNumVars; i++) { s.newVar(); }

        AND_GIVEN("a typical extension variable definition") {
            std::vector<CRef> cs;
            Lit a = mkLit(1), b = mkLit(2), x = mkLit(3);
            setVec(ps, {~x, a, b}); cs.push_back(s.ca.alloc(ps)); s.attachClause(cs[cs.size() - 1]);
            setVec(ps, { x,~a   }); cs.push_back(s.ca.alloc(ps)); s.attachClause(cs[cs.size() - 1]);
            setVec(ps, { x,   ~b}); cs.push_back(s.ca.alloc(ps)); s.attachClause(cs[cs.size() - 1]);

            WHEN("the basis literals are falsified") {
                ser.set_value(var(a), l_False, 0);
                ser.set_value(var(b), l_False, 1);

                cr = ser.findAssertingClause(i_undef, i_max, x, cs);
                REQUIRE(cr == cs[0]);
                CHECK(i_undef == 0 );
                CHECK(i_max   == 2 );
            }

            WHEN("the basis literals are satisfied") {
                ser.set_value(var(a), l_True, 0);
                ser.set_value(var(b), l_True, 1);

                cr = ser.findAssertingClause(i_undef, i_max, ~x, cs);
                REQUIRE((cr == cs[1] || cr == cs[2]));
                CHECK(i_undef == 0 );
                CHECK(i_max   == 1 );
            }
        }

        AND_GIVEN("any set of clauses containing exactly one asserting clause") {
            std::vector<CRef> cs;

            // Non-asserting clauses
            ps.clear(); for (int i = 0; i < ser.originalNumVars; i++) ps.push(mkLit(i, i % 2 == 0)); 
            cs.push_back(s.ca.alloc(ps)); s.attachClause(cs[cs.size() - 1]);
            ps.clear(); for (int i = 0; i < ser.originalNumVars; i++) ps.push(mkLit(i, true)); 
            cs.push_back(s.ca.alloc(ps)); s.attachClause(cs[cs.size() - 1]);

            // Asserting clause
            int asserting_cr = cs.size();
            ps.clear(); for (int i = 0; i < ser.originalNumVars; i++) ps.push(mkLit(i, false)); 
            cs.push_back(s.ca.alloc(ps)); s.attachClause(cs[cs.size() - 1]);

            // Non-asserting clauses
            ps.clear(); for (int i = 0; i < ser.originalNumVars; i++) ps.push(mkLit(i, i % 2 == 1)); 
            cs.push_back(s.ca.alloc(ps)); s.attachClause(cs[cs.size() - 1]);
            ps.clear(); for (int i = 0; i < ser.originalNumVars; i++) ps.push(mkLit(i, i % 3 == 0)); 
            cs.push_back(s.ca.alloc(ps)); s.attachClause(cs[cs.size() - 1]);

            WHEN("the asserting literal is at the beginning") {
                int expect_i_undef = 0, expect_i_max = 4;
                setVariables(ser, expect_i_undef, expect_i_max, ser.originalNumVars);

                Lit x = s.ca[cs[asserting_cr]][expect_i_undef];
                cr = ser.findAssertingClause(i_undef, i_max, ~x, cs);
                
                THEN("the correct clause and indices are found") {
                    REQUIRE(cr == cs[asserting_cr]);
                    CHECK(i_undef == expect_i_undef);
                    CHECK(i_max   == expect_i_max);
                }
            }

            WHEN("the asserting literal is at the end") {
                int expect_i_undef = ser.originalNumVars - 1, expect_i_max = 4;
                setVariables(ser, expect_i_undef, expect_i_max, ser.originalNumVars);
                Lit x = s.ca[cs[asserting_cr]][expect_i_undef];
                cr = ser.findAssertingClause(i_undef, i_max, ~x, cs);
                
                THEN("the correct clause and indices are found") {
                    REQUIRE(cr == cs[asserting_cr]);
                    CHECK(i_undef == expect_i_undef);
                    CHECK(i_max   == expect_i_max);
                }
            }

            WHEN("the highest-level literal is at the beginning") {
                int expect_i_undef = 4, expect_i_max = 0;
                setVariables(ser, expect_i_undef, expect_i_max, ser.originalNumVars);
                Lit x = s.ca[cs[asserting_cr]][expect_i_undef];
                cr = ser.findAssertingClause(i_undef, i_max, ~x, cs);
                
                THEN("the correct clause and indices are found") {
                    REQUIRE(cr == cs[asserting_cr]);
                    CHECK(i_undef == expect_i_undef);
                    CHECK(i_max   == expect_i_max);
                }
            }

            WHEN("the highest-level literal is at the end") {
                int expect_i_undef = 4, expect_i_max = ser.originalNumVars - 1;
                setVariables(ser, expect_i_undef, expect_i_max, ser.originalNumVars);
                Lit x = s.ca[cs[asserting_cr]][expect_i_undef];
                cr = ser.findAssertingClause(i_undef, i_max, ~x, cs);
                
                THEN("the correct clause and indices are found") {
                    REQUIRE(cr == cs[asserting_cr]);
                    CHECK(i_undef == expect_i_undef);
                    CHECK(i_max   == expect_i_max);
                }
            }

            WHEN("the asserting literal occurs before the highest-level literal") {
                int expect_i_undef = 3, expect_i_max = 5;
                setVariables(ser, expect_i_undef, expect_i_max, ser.originalNumVars);
                Lit x = s.ca[cs[asserting_cr]][expect_i_undef];
                cr = ser.findAssertingClause(i_undef, i_max, ~x, cs);
                
                THEN("the correct clause and indices are found") {
                    REQUIRE(cr == cs[asserting_cr]);
                    CHECK(i_undef == expect_i_undef);
                    CHECK(i_max   == expect_i_max);
                }
            }

            WHEN("the asserting literal occurs after the highest-level literal") {
                int expect_i_undef = 5, expect_i_max = 3;
                setVariables(ser, expect_i_undef, expect_i_max, ser.originalNumVars);
                Lit x = s.ca[cs[asserting_cr]][expect_i_undef];
                cr = ser.findAssertingClause(i_undef, i_max, ~x, cs);
                
                THEN("the correct clause and indices are found") {
                    REQUIRE(cr == cs[asserting_cr]);
                    CHECK(i_undef == expect_i_undef);
                    CHECK(i_max   == expect_i_max);
                }
            }
        }
    }
}

// SolverER::deleteExtVars
// SCENARIO("Deleting extension variables", "[SolverER]") {}
}