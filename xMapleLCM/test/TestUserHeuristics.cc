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