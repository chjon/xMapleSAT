/**********************************************************************[RestartHeuristicManager.cc]
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

#include "core/RestartHeuristicManager.h"
#include "core/Solver.h"

using namespace Minisat;

///////////////////////////////////////////////////////////////////////////////////////////////////
// OPTIONS

static const char* _cat = "CORE";

static IntOption    opt_restart_first(_cat, "rfirst",      "The base restart interval", 100, IntRange(1, INT32_MAX));
static DoubleOption opt_restart_inc  (_cat, "rinc",        "Restart interval increase factor", 2, DoubleRange(1, false, HUGE_VAL, false));

///////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS

RestartHeuristicManager::RestartHeuristicManager(Solver& s)
    : VSIDS(s.branchingHeuristicManager.VSIDS)
    , cached(false)
    , global_lbd_sum(0)
    , lbd_queue(50)
    , conflictBudget(0)

    , restart_first(opt_restart_first)
    , restart_inc(opt_restart_inc)

    , curr_restarts(0)
    , conflicts_VSIDS(0)
{}
