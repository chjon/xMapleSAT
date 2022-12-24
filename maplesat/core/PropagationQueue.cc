#include "core/PropagationQueue.h"
#include "core/Solver.h"

using namespace Minisat;

PropagationQueue::PropagationQueue(Solver* s)
#if BCP_PRIORITY_MODE == BCP_PRIORITY_IMMEDIATE
    : qhead(0)
    , queue(s->assignmentTrail.trail) 

#elif BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
    : qhead(0)
    , queue(s->assignmentTrail.trail) 
    , order_heap(LitOrderLt<double>(s->branchingHeuristicManager.getActivityVSIDS()))

#elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
    : order_heap(LitOrderLt<double>(s->branchingHeuristicManager.getActivityVSIDS()))

#endif
    , variableDatabase(s->variableDatabase)
    , assignmentTrail(s->assignmentTrail)
{}