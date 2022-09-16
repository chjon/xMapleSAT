#include "catch.hpp"
#include <core/SolverER.h>
#include <test/Util.h>

namespace Minisat {
SCENARIO("Generating definitions for extension variables", "[UserHeuristics]") {
    GIVEN("clauses with sufficiently many variables") {
        Solver s;
        SolverER& ser = *(s.ser);
        std::vector<ExtDef>& output = ser.m_extVarDefBuffer;
        std::vector<CRef>& selectedClauses = ser.m_selectedClauses;

        // Set up some input clauses
        vec<Lit> ps;
        for (int i = 0; i < 10; i++) s.newVar();
        setLitVec(ps, { 1, 2, 3, 4}); s.addClause_(ps);
        setLitVec(ps, {-1, 2,   -4}); s.addClause_(ps);
        setLitVec(ps, {-1,-2,-3,-4}); s.addClause_(ps);
        setLitVec(ps, {    2, 3,-4}); s.addClause_(ps);
        setLitVec(ps, { 1,    3,-4}); s.addClause_(ps);
        s.ser->originalNumVars = s.nVars();

        // Select all the clauses
        for (int i = 0; i < s.clauses.size(); i++) {
            selectedClauses.push_back(s.clauses[i]);
        }

        WHEN("the requested number of definitions is less than or equal to the number of possible definitions") {
            const unsigned int maxNumNewVars = 5;
            AND_WHEN("using user_extDefHeuristic_random") {
                ser.user_extDefHeuristic_random(output, selectedClauses, maxNumNewVars);

                THEN("the requested number of variables should be generated") {
                    CHECK(output.size() == maxNumNewVars);
                }

                THEN("no duplicates should be generated") {
                    CHECK_THAT(output, extDefUnique());
                }
            }

            AND_WHEN("using user_extDefHeuristic_subexpression") {
                ser.user_extDefHeuristic_subexpression(output, selectedClauses, maxNumNewVars);

                THEN("the requested number of variables should be generated") {
                    CHECK(output.size() == maxNumNewVars);
                }

                THEN("no duplicates should be generated") {
                    CHECK_THAT(output, extDefUnique());
                }
            }
        }

        WHEN("the requested number of definitions exceeds the number of possible definitions") {
            const unsigned int maxNumNewVars = 50;
            AND_WHEN("using user_extDefHeuristic_random") {
                ser.user_extDefHeuristic_random(output, selectedClauses, maxNumNewVars);

                THEN("no duplicates should be generated") {
                    CHECK_THAT(output, extDefUnique());
                }
            }

            AND_WHEN("using user_extDefHeuristic_subexpression") {
                ser.user_extDefHeuristic_subexpression(output, selectedClauses, maxNumNewVars);

                THEN("no duplicates should be generated") {
                    CHECK_THAT(output, extDefUnique());
                }
            }
        }
    }
}
}