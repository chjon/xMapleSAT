/******************************************************************************************[ExtDefMap.h]
xMaple* -- Copyright (c) 2022, Jonathan Chung, Vijay Ganesh, Sam Buss

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

#include "catch.hpp"
#include "core/SolverTypes.h"
#include "er/ERSolver.h"
#include "mtl/Vec.h"
#include "test/Util.h"

namespace Minisat {

SCENARIO("Introducing extension variables", "[ERManager]") {
    GIVEN("basis variables") {
        ERSolver s;
        AssignmentTrail& at = s.assignmentTrail;
        ERManager& erm = s.erManager;
        std::tr1::unordered_map<Var, std::vector<CRef> > db;
        std::vector< std::vector<Lit> > additional;
        vec<Lit> clause, expect;

        // Set up variables for testing
        erm.originalNumVars = 10;
        for (int i = 0; i < erm.originalNumVars; i++) { s.newVar(); }

        WHEN("introducing a new definition") {
            Lit x = mkLit(10), a = mkLit(0), b = mkLit(1);
            additional.push_back(std::vector<Lit>({x, a, b, mkLit(2)}));
            erm.m_extVarDefBuffer.push_back(ExtDef{ x, a, b, additional });
            erm.introduceExtVars(db);

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
                REQUIRE(at.nVars() == erm.originalNumVars + 1);

                // Test whether the extension variable definition was stored in the extension definition map
                REQUIRE(erm.isCurrentExtVar(var(x)));

                // Check whether the buffer was cleared
                REQUIRE(erm.m_extVarDefBuffer.size() == 0);
            }
        }
    }
}

SCENARIO("Testing for valid definition pairs", "[ERManager]") {
    ERSolver s;
    AssignmentTrail& at = s.assignmentTrail;
    ERManager& erm = s.erManager;
    std::tr1::unordered_map<Var, std::vector<CRef> > db;

    // Set up variables for testing
    erm.originalNumVars = 10;
    for (int i = 0; i < erm.originalNumVars; i++) { s.newVar(); }

    GIVEN("no pre-existing pairs") {
        std::tr1::unordered_set< std::pair<Lit, Lit> > generatedPairs;

        WHEN("the input pair is two of the same variable") {
            Lit a = mkLit(1);

            THEN("the pair is not valid") {
                CHECK_FALSE(erm.isValidDefPair( a,  a, generatedPairs));
                CHECK_FALSE(erm.isValidDefPair( a, ~a, generatedPairs));
                CHECK_FALSE(erm.isValidDefPair(~a,  a, generatedPairs));
                CHECK_FALSE(erm.isValidDefPair(~a, ~a, generatedPairs));
            }
        }

        WHEN("the input pair consists of different variables") {
            Lit a = mkLit(1, false), b = mkLit(2, false);

            THEN("the pair is valid") {
                CHECK(erm.isValidDefPair( a,  b, generatedPairs));
                CHECK(erm.isValidDefPair(~a,  b, generatedPairs));
                CHECK(erm.isValidDefPair(~a, ~b, generatedPairs));
                CHECK(erm.isValidDefPair( a, ~b, generatedPairs));
            }
        }

        WHEN("checking input pairs with different assignment levels") {
            at.assign(mkLit(1, false)); at.assign(mkLit(2, true)); at.newDecisionLevel();
            at.assign(mkLit(3, false)); at.assign(mkLit(4, true));

            THEN("pairs containing a literal set at level 0 should be rejected") {
                CHECK_FALSE(erm.isValidDefPair(mkLit(1), mkLit(2), generatedPairs));
                CHECK_FALSE(erm.isValidDefPair(mkLit(1), mkLit(3), generatedPairs));
                CHECK_FALSE(erm.isValidDefPair(mkLit(1), mkLit(4), generatedPairs));
                CHECK_FALSE(erm.isValidDefPair(mkLit(3), mkLit(2), generatedPairs));
                CHECK_FALSE(erm.isValidDefPair(mkLit(4), mkLit(2), generatedPairs));
                CHECK      (erm.isValidDefPair(mkLit(3), mkLit(4), generatedPairs));
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
                CHECK(erm.isValidDefPair(a, b, generatedPairs));
                CHECK(erm.isValidDefPair(b, a, generatedPairs));
            }
        }

        WHEN("the input pair has single overlap with existing pairs") {
            Lit a = mkLit(1), b = mkLit(6);
            
            THEN("the pair should be accepted") {
                CHECK(erm.isValidDefPair(a, b, generatedPairs));
                CHECK(erm.isValidDefPair(b, a, generatedPairs));
            }
        }

        WHEN("the input pair has double overlap with existing pairs") {
            Lit a = mkLit(1), b = mkLit(3);
            
            THEN("the pair should be accepted") {
                CHECK(erm.isValidDefPair(a, b, generatedPairs));
                CHECK(erm.isValidDefPair(b, a, generatedPairs));
            }
        }

        WHEN("the input pair has already been added") {
            Lit a = mkLit(1), b = mkLit(2);
            
            THEN("the pair should be rejected") {
                CHECK_FALSE(erm.isValidDefPair(a, b, generatedPairs));
                CHECK_FALSE(erm.isValidDefPair(b, a, generatedPairs));
            }
        }
    }
}

SCENARIO("Choosing extension variables to delete", "[ERManager]") {

    GIVEN("Extension variables") {
        using DeletionPredicate = std::function<bool(Var)>;
        ERSolver s;
        ERManager& erm = s.erManager;
        std::vector< std::vector<Lit> > additional;
        std::tr1::unordered_set<Var> varsToDelete;
        vec<Var> actual, expect;

        // Set up variables for testing
        erm.originalNumVars = 10;
        for (int i = 0; i < erm.originalNumVars; i++) { s.newVar(); }

        erm.m_extVarDefBuffer.push_back(ExtDef{ mkLit(10), mkLit(1), mkLit(2), additional });
        erm.m_extVarDefBuffer.push_back(ExtDef{ mkLit(11), mkLit(3), mkLit(4), additional });
        erm.m_extVarDefBuffer.push_back(ExtDef{ mkLit(12), mkLit(5), mkLit(6), additional });
        erm.m_extVarDefBuffer.push_back(ExtDef{ mkLit(13), mkLit(1), mkLit(4), additional });
        erm.introduceExtVars(erm.extDefs);

        WHEN("there are no restrictions") {
            DeletionPredicate deletionPredicate = [](Var v){ return true; };
            erm.getExtVarsToDelete(varsToDelete, deletionPredicate);

            THEN("all extension variables are chosen for deletion") {
                setVec(actual, varsToDelete);
                setVec(expect, {10, 11, 12, 13});
                REQUIRE_THAT(actual, vecEqualUnordered(expect));
            }
        }

        WHEN("some extension variables are basis literals") {
            DeletionPredicate deletionPredicate = [](Var v){ return true; };

            erm.m_extVarDefBuffer.push_back(ExtDef{ mkLit(14), mkLit( 3), mkLit(10), additional });
            erm.m_extVarDefBuffer.push_back(ExtDef{ mkLit(15), mkLit(11), mkLit(12), additional });
            erm.introduceExtVars(erm.extDefs);
            erm.getExtVarsToDelete(varsToDelete, deletionPredicate);

            THEN("only non-basis literal extension variables are chosen for deletion") {
                setVec(actual, varsToDelete);
                setVec(expect, {13, 14, 15});
                REQUIRE_THAT(actual, vecEqualUnordered(expect));
            }
        }

        WHEN("the deletion predicate rejects every variable") {
            DeletionPredicate deletionPredicate = [](Var v){ return false; };
            erm.getExtVarsToDelete(varsToDelete, deletionPredicate);

            THEN("no extension variables are chosen for deletion") {
                setVec(actual, varsToDelete);
                setVec(expect, {});
                REQUIRE_THAT(actual, vecEqualUnordered(expect));
            }
        }

        WHEN("the deletion predicate rejects some variables") {
            DeletionPredicate deletionPredicate = [](Var v){ return v % 2 == 0; };
            erm.getExtVarsToDelete(varsToDelete, deletionPredicate);

            THEN("none of the rejected variables are selected") {
                setVec(actual, varsToDelete);
                setVec(expect, {10, 12});
                REQUIRE_THAT(actual, vecEqualUnordered(expect));
            }
        }
    }
}

SCENARIO("Find asserting extension definition clause", "[ERManager]") {
    GIVEN("Variables") {
        ERSolver s;
        AssignmentTrail& at = s.assignmentTrail;
        ERManager& erm = s.erManager;
        UnitPropagator& up = s.unitPropagator;
        std::tr1::unordered_set<Var> varsToDelete;
        vec<Lit> ps, actual, expect;
        int i_undef = -1, i_max = -1;
        CRef cr;

        // Set up variables for testing
        erm.originalNumVars = 10;
        for (int i = 0; i < erm.originalNumVars; i++) { s.newVar(); }

        AND_GIVEN("a typical extension variable definition") {
            std::vector<CRef> cs;
            Lit a = mkLit(1), b = mkLit(2), x = mkLit(3);
            setVec(ps, {~x, a, b}); cs.push_back(s.ca.alloc(ps)); up.attachClause(cs[cs.size() - 1]);
            setVec(ps, { x,~a   }); cs.push_back(s.ca.alloc(ps)); up.attachClause(cs[cs.size() - 1]);
            setVec(ps, { x,   ~b}); cs.push_back(s.ca.alloc(ps)); up.attachClause(cs[cs.size() - 1]);

            WHEN("the basis literals are falsified") {
                at.cancelUntil(0);
                at.newDecisionLevel(); at.assign(~a); 
                at.newDecisionLevel(); at.assign(~b);

                cr = erm.findAssertingClause(i_undef, i_max, x, cs);
                REQUIRE(cr == cs[0]);
                CHECK(i_undef == 0 );
                CHECK(i_max   == 2 );
            }

            WHEN("the basis literals are satisfied") {
                at.cancelUntil(0);
                at.newDecisionLevel(); at.assign(a);
                at.newDecisionLevel(); at.assign(b);

                cr = erm.findAssertingClause(i_undef, i_max, ~x, cs);
                REQUIRE((cr == cs[1] || cr == cs[2]));
                CHECK(i_undef == 0 );
                CHECK(i_max   == 1 );
            }
        }

        AND_GIVEN("any set of clauses containing exactly one asserting clause") {
            std::vector<CRef> cs;

            // Non-asserting clauses
            ps.clear(); for (int i = 0; i < erm.originalNumVars; i++) ps.push(mkLit(i, i % 2 == 0)); 
            cs.push_back(s.ca.alloc(ps)); up.attachClause(cs[cs.size() - 1]);
            ps.clear(); for (int i = 0; i < erm.originalNumVars; i++) ps.push(mkLit(i, true)); 
            cs.push_back(s.ca.alloc(ps)); up.attachClause(cs[cs.size() - 1]);

            // Asserting clause
            int asserting_cr = cs.size();
            ps.clear(); for (int i = 0; i < erm.originalNumVars; i++) ps.push(mkLit(i, false)); 
            cs.push_back(s.ca.alloc(ps)); up.attachClause(cs[cs.size() - 1]);

            // Non-asserting clauses
            ps.clear(); for (int i = 0; i < erm.originalNumVars; i++) ps.push(mkLit(i, i % 2 == 1)); 
            cs.push_back(s.ca.alloc(ps)); up.attachClause(cs[cs.size() - 1]);
            ps.clear(); for (int i = 0; i < erm.originalNumVars; i++) ps.push(mkLit(i, i % 3 == 0)); 
            cs.push_back(s.ca.alloc(ps)); up.attachClause(cs[cs.size() - 1]);

            WHEN("the asserting literal is at the beginning") {
                int expect_i_undef = 0, expect_i_max = 4;
                setVariables(at, expect_i_undef, expect_i_max, erm.originalNumVars);

                Lit x = s.ca[cs[asserting_cr]][expect_i_undef];
                cr = erm.findAssertingClause(i_undef, i_max, ~x, cs);
                
                THEN("the correct clause and indices are found") {
                    REQUIRE(cr == cs[asserting_cr]);
                    CHECK(i_undef == expect_i_undef);
                    CHECK(i_max   == expect_i_max);
                }
            }

            WHEN("the asserting literal is at the end") {
                int expect_i_undef = erm.originalNumVars - 1, expect_i_max = 4;
                setVariables(at, expect_i_undef, expect_i_max, erm.originalNumVars);
                Lit x = s.ca[cs[asserting_cr]][expect_i_undef];
                cr = erm.findAssertingClause(i_undef, i_max, ~x, cs);
                
                THEN("the correct clause and indices are found") {
                    REQUIRE(cr == cs[asserting_cr]);
                    CHECK(i_undef == expect_i_undef);
                    CHECK(i_max   == expect_i_max);
                }
            }

            WHEN("the highest-level literal is at the beginning") {
                int expect_i_undef = 4, expect_i_max = 0;
                setVariables(at, expect_i_undef, expect_i_max, erm.originalNumVars);
                Lit x = s.ca[cs[asserting_cr]][expect_i_undef];
                cr = erm.findAssertingClause(i_undef, i_max, ~x, cs);
                
                THEN("the correct clause and indices are found") {
                    REQUIRE(cr == cs[asserting_cr]);
                    CHECK(i_undef == expect_i_undef);
                    CHECK(i_max   == expect_i_max);
                }
            }

            WHEN("the highest-level literal is at the end") {
                int expect_i_undef = 4, expect_i_max = erm.originalNumVars - 1;
                setVariables(at, expect_i_undef, expect_i_max, erm.originalNumVars);
                Lit x = s.ca[cs[asserting_cr]][expect_i_undef];
                cr = erm.findAssertingClause(i_undef, i_max, ~x, cs);
                
                THEN("the correct clause and indices are found") {
                    REQUIRE(cr == cs[asserting_cr]);
                    CHECK(i_undef == expect_i_undef);
                    CHECK(i_max   == expect_i_max);
                }
            }

            WHEN("the asserting literal occurs before the highest-level literal") {
                int expect_i_undef = 3, expect_i_max = 5;
                setVariables(at, expect_i_undef, expect_i_max, erm.originalNumVars);
                Lit x = s.ca[cs[asserting_cr]][expect_i_undef];
                cr = erm.findAssertingClause(i_undef, i_max, ~x, cs);
                
                THEN("the correct clause and indices are found") {
                    REQUIRE(cr == cs[asserting_cr]);
                    CHECK(i_undef == expect_i_undef);
                    CHECK(i_max   == expect_i_max);
                }
            }

            WHEN("the asserting literal occurs after the highest-level literal") {
                int expect_i_undef = 5, expect_i_max = 3;
                setVariables(at, expect_i_undef, expect_i_max, erm.originalNumVars);
                Lit x = s.ca[cs[asserting_cr]][expect_i_undef];
                cr = erm.findAssertingClause(i_undef, i_max, ~x, cs);
                
                THEN("the correct clause and indices are found") {
                    REQUIRE(cr == cs[asserting_cr]);
                    CHECK(i_undef == expect_i_undef);
                    CHECK(i_max   == expect_i_max);
                }
            }
        }
    }
}

// ERManager::deleteExtVars
// SCENARIO("Deleting extension variables", "[ERManager]") {}
}