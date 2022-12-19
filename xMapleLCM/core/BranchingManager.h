#ifndef Minisat_BranchingManager_h
#define Minisat_BranchingManager_h

#include "core/SolverTypes.h"

namespace Minisat {
    class BranchingManager {
        void     insertVarOrder   (Var x); // Insert a variable in the decision order priority queue.
        Lit      pickBranchLit    ();      // Return the next decision variable.
        void     unsetVar         ();      // Unassign a variable and put it back in the priority queue.
    };
} // namespace Minisat

#endif