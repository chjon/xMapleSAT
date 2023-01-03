/*****************************************************************************[PropagationQueue.cc]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010, Niklas Sorensson

MapleSAT_Refactor, based on MapleSAT -- Copyright (c) 2022, Jonathan Chung, Vijay Ganesh, Sam Buss

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

#include "core/PropagationQueue.h"
#include "core/Solver.h"

using namespace Minisat;

PropagationQueue::PropagationQueue(Solver& s)
#if BCP_PRIORITY_MODE == BCP_PRIORITY_IMMEDIATE
    : qhead(0)
    , queue(s.assignmentTrail.getTrail()) 

#elif BCP_PRIORITY_MODE == BCP_PRIORITY_DELAYED
    : qhead(0)
    , queue(s.assignmentTrail.getTrail()) 
    , order_heap(LitOrderLt<double>(s.branchingHeuristicManager.getActivityVSIDS()))

#elif BCP_PRIORITY_MODE == BCP_PRIORITY_OUT_OF_ORDER
    : order_heap(LitOrderLt<double>(s.branchingHeuristicManager.getActivityVSIDS()))

#endif
    , variableDatabase(s.variableDatabase)
    , assignmentTrail(s.assignmentTrail)
{}