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
#include <test/Util.h>
#include <core/Solver.h>
#include <core/SolverTypes.h>
#include <mtl/Vec.h>

namespace Minisat {

SCENARIO("Enforce watcher invariant", "[PropagationComponent]") {
    GIVEN("A clause") {
        Solver s;
        PropagationComponent& pc = s.propagationComponent;
        std::tr1::unordered_set<Lit> varsToDelete;
        vec<Lit> ps, actual, expect;

        // Set up variables for testing
        int originalNumVars = 5;
        for (int i = 0; i < originalNumVars; i++) { s.newVar(); }

        // Set up clause
        ps.clear(); for (int i = 0; i < originalNumVars; i++) ps.push(mkLit(i)); 
        CRef cr = s.ca.alloc(ps); s.attachClause(cr);
        Clause& c = s.ca[cr];

        // Ensure watchers are set up correctly
        clause2Vec(actual, c);
        REQUIRE_THAT(actual, watchersCorrect(pc.watches, cr));

        WHEN("the asserting literal is at index 0 and the highest-level literal is at index 1") {
            int i_undef = 0, i_max = 1;
            setVariables(pc, i_undef, i_max, originalNumVars);
            pc.enforceWatcherInvariant(cr, i_undef, i_max);

            THEN("the clause does not change and the watchers are correct") {
                setVec(expect, {mkLit(0), mkLit(1), mkLit(2), mkLit(3), mkLit(4)});
                clause2Vec(actual, c);
                REQUIRE_THAT(actual, vecEqual(expect));
                REQUIRE_THAT(actual, watchersCorrect(pc.watches, cr));
            }
        }

        WHEN("the asserting literal is at index 1 and the highest-level literal is at index 0") {
            int i_undef = 1, i_max = 0;
            setVariables(pc, i_undef, i_max, originalNumVars);
            pc.enforceWatcherInvariant(cr, i_undef, i_max);

            THEN("the first two literals in the clause are swapped and the watchers are correct") {
                setVec(expect, {mkLit(1), mkLit(0), mkLit(2), mkLit(3), mkLit(4)});
                clause2Vec(actual, c);
                REQUIRE_THAT(actual, vecEqual(expect));
                REQUIRE_THAT(actual, watchersCorrect(pc.watches, cr));
            }
        }

        WHEN("the asserting literal is at index 0 and the highest-level literal is at index > 1") {
            int i_undef = 0, i_max = 4;
            setVariables(pc, i_undef, i_max, originalNumVars);
            pc.enforceWatcherInvariant(cr, i_undef, i_max);

            THEN("the highest-level literal is moved to index 1 and the watchers are correct") {
                setVec(expect, {mkLit(0), mkLit(4), mkLit(2), mkLit(3), mkLit(1)});
                clause2Vec(actual, c);
                REQUIRE_THAT(actual, vecEqual(expect));
                REQUIRE_THAT(actual, watchersCorrect(pc.watches, cr));
            }
        }

        WHEN("the highest-level literal is at index 1 and the asserting literal is at index > 1") {
            int i_undef = 4, i_max = 1;
            setVariables(pc, i_undef, i_max, originalNumVars);
            pc.enforceWatcherInvariant(cr, i_undef, i_max);

            THEN("the asserting literal is moved to index 0 and the watchers are correct") {
                setVec(expect, {mkLit(4), mkLit(1), mkLit(2), mkLit(3), mkLit(0)});
                clause2Vec(actual, c);
                REQUIRE_THAT(actual, vecEqual(expect));
                REQUIRE_THAT(actual, watchersCorrect(pc.watches, cr));
            }
        }

        WHEN("the highest-level literal is at index 0 and the asserting literal is at index > 1") {
            int i_undef = 4, i_max = 0;
            setVariables(pc, i_undef, i_max, originalNumVars);
            pc.enforceWatcherInvariant(cr, i_undef, i_max);

            THEN("the asserting literal is moved to index 0, the highest-level literal is moved to index 1, and the watchers are correct") {
                setVec(expect, {mkLit(4), mkLit(0), mkLit(2), mkLit(3), mkLit(1)});
                clause2Vec(actual, c);
                REQUIRE_THAT(actual, vecEqual(expect));
                REQUIRE_THAT(actual, watchersCorrect(pc.watches, cr));
            }
        }

        WHEN("the asserting literal is at index 1 and the highest-level literal is at index > 1") {
            int i_undef = 1, i_max = 4;
            setVariables(pc, i_undef, i_max, originalNumVars);
            pc.enforceWatcherInvariant(cr, i_undef, i_max);

            THEN("the asserting literal is moved to index 0, the highest-level literal is moved to index 1, and the watchers are correct") {
                setVec(expect, {mkLit(1), mkLit(4), mkLit(2), mkLit(3), mkLit(0)});
                clause2Vec(actual, c);
                REQUIRE_THAT(actual, vecEqual(expect));
                REQUIRE_THAT(actual, watchersCorrect(pc.watches, cr));
            }
        }

        WHEN("the asserting literal is at index > 1 and the highest-level literal is at index > 1") {
            int i_undef = 4, i_max = 3;
            setVariables(pc, i_undef, i_max, originalNumVars);
            pc.enforceWatcherInvariant(cr, i_undef, i_max);

            THEN("the asserting literal is moved to index 0, the highest-level literal is moved to index 1, and the watchers are correct") {
                setVec(expect, {mkLit(4), mkLit(3), mkLit(2), mkLit(1), mkLit(0)});
                clause2Vec(actual, c);
                REQUIRE_THAT(actual, vecEqual(expect));
                REQUIRE_THAT(actual, watchersCorrect(pc.watches,cr));
            }
        }
    }
}

// SolverER::deleteExtVars
// SCENARIO("Deleting extension variables", "[SolverER]") {}
}