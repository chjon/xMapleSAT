# Adding the Extended Resolution Framework to Maplesat Variants
## Changes to Solver.h
- additional #includes:
    - `#include "core/SolverERTypes.h"`
- forward declarations:
    - add `class SolverER;` before the main `Solver` class
- class modifications
    - add `SolverER* ser;` as a public class member
    - `friend class SolverER;`

## Changes to Solver.cc
- additional #includes:
    - `#include "core/SolverER.h"`
- constructor/destructor:
    - add `ser = new SolverER(this);` to the constructor
    - add `delete ser;` to the destructor
- changes to `Solver::simplify`:
    - call `ser->removeSatisfied();` similarly to `removeSatisfied`

- changes to `Solver::search`:
    - call `ser->generateDefinitions();`, then `ser->introduceExtVars(ser->extDefs);` before the main search loop
    - call `ser->substitute(learnt_clause, ser->user_extSubPredicate);` after backtracking and before adding the clause to the learnt clause database
    - call `ser->deleteExtVars(ser->user_extDelPredicate);` immediately before learnt clause deletion (`Solver::reduceDB`)

    - statistics: after picking a branch literal, add `if (ser->isExtVar(var(next))) ser->branchOnExt++;`

- changes to `Solver::solve_`:
    - add `ser->originalNumVars = nVars();` before the loop that calls `search`

- changes to `Solver::relocAll`:
    - call `ser->relocAll(to);` similarly to the original and learnt clause databases