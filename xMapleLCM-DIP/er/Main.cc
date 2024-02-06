/*****************************************************************************************[Main.cc]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010, Niklas Sorensson

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

#include <errno.h>

#include <signal.h>
#include <zlib.h>

#include "utils/System.h"
#include "utils/ParseUtils.h"
#include "utils/Options.h"
#include "core/Dimacs.h"
#include "er/ERSolver.h"

using namespace Minisat;

//=================================================================================================


void printStats(ERSolver& solver) {
    double cpu_time = cpuTime();
    double mem_used = memUsedPeak();
    const ERManager& erm = solver.erManager;

    // Standard stats
    printf("restarts                  : %"    PRIu64 "\n", solver.starts);
    printf("conflicts                 : %-12" PRIu64 "   (%.0f /sec)\n", solver.conflicts   , solver.conflicts   /cpu_time);
    printf("decisions                 : %-12" PRIu64 "   (%4.2f %% random) (%.0f /sec)\n", solver.branchingHeuristicManager.decisions, 0.f, solver.branchingHeuristicManager.decisions   /cpu_time);
    printf("propagations              : %-12" PRIu64 "   (%.0f /sec)\n", solver.unitPropagator.propagations, solver.unitPropagator.propagations/cpu_time);
    printf("conflict literals         : %-12" PRIu64 "   (%4.2f %% deleted)\n", solver.conflictAnalyzer.tot_literals, (solver.conflictAnalyzer.max_literals - solver.conflictAnalyzer.tot_literals)*100 / (double)solver.conflictAnalyzer.max_literals);
    printf("\n");
    
    // Extended resolution stats
    printf("conflicts w/ DIP found    : %-12" PRIu64 "   (%4.2f %% of conflicts)\n", solver.conflictAnalyzer.conflicts_with_dip, double(solver.conflictAnalyzer.conflicts_with_dip)/solver.conflicts*100);
    printf("conflicts w/ dangerous DIP: %-12" PRIu64 "   (%4.2f %% of conflicts)\n", solver.conflictAnalyzer.conflicts_with_dangerous_dip, double(solver.conflictAnalyzer.conflicts_with_dangerous_dip)/solver.conflicts*100);
    printf("conflicts w/ DIP-learning : %-12" PRIu64 "   (%4.2f %% of conflicts)\n", solver.dip_conflicts, double(solver.dip_conflicts)/solver.conflicts*100);
    printf("decisions on ext vars     : %-12" PRIu64 "   (%4.2f %% of decisions)\n", erm.branchOnExt, double(erm.branchOnExt)/solver.branchingHeuristicManager.decisions*100);

    printf("total ext vars            : %-12" PRIu64 "\n", erm.total_ext_vars);
    printf("tried delete ext vars     : %-12" PRIu64 "\n", erm.tried_del_ext_vars);
    printf("deleted ext vars          : %-12" PRIu64 "\n", erm.deleted_ext_vars);
    printf("max ext vars              : %-12" PRIu64 "\n", erm.max_ext_vars);
    // printf("conflict ext clauses   : %-12" PRIu64 "   (%.0f /sec)\n", erm.conflict_extclauses, erm.conflict_extclauses / cpu_time);
    // printf("learnt ext clauses     : %-12" PRIu64 "   (%.0f /sec)\n", erm.learnt_extclauses, erm.learnt_extclauses / cpu_time);

    printf("\n");
    // Resource usage stats
    if (mem_used != 0) printf("Memory used               : %.2f MB\n", mem_used);
    printf("CPU time                  : %g s\n", cpu_time);

    // ER overhead stats
    printf("DIP computation time      : %g s (%4.2f %% of total time)\n", solver.conflictAnalyzer.time_DIP,solver.conflictAnalyzer.time_DIP/cpu_time*100);
    printf("ER_sel time               : %g s\n", erm.extTimerRead(0));
    printf("ER_add time               : %g s\n", erm.extTimerRead(1));
    printf("ER_delC time              : %g s\n", erm.extTimerRead(2));
    printf("ER_delV time              : %g s\n", erm.extTimerRead(3));
    printf("ER_sub time               : %g s\n", erm.extTimerRead(4));
    printf("ER_stat time              : %g s\n", erm.extTimerRead(5));
}


static ERSolver* solver;
// Terminate by notifying the solver and back out gracefully. This is mainly to have a test-case
// for this feature of the Solver as it may take longer than an immediate call to '_exit()'.
static void SIGINT_interrupt(int signum) { solver->interrupt(); }

// Note that '_exit()' rather than 'exit()' has to be used. The reason is that 'exit()' calls
// destructors and may cause deadlocks if a malloc/free function happens to be running (these
// functions are guarded by locks for multithreaded use).
static void SIGINT_exit(int signum) {
    printf("\n"); printf("*** INTERRUPTED ***\n");
    if (solver->verbosity > 0){
        printStats(*solver);
        printf("\n"); printf("*** INTERRUPTED ***\n"); }
    _exit(1); }


//=================================================================================================
// Main:


int main(int argc, char** argv) {
    try {
        setUsageHelp("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
        // printf("This is MiniSat 2.0 beta\n");
        
#if defined(__linux__) && defined(_FPU_EXTENDED) && defined(_FPU_DOUBLE) && defined(_FPU_GETCW)
        fpu_control_t oldcw, newcw;
        _FPU_GETCW(oldcw); newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(newcw);
        printf("WARNING: for repeatability, setting FPU to use double precision\n");
#endif
        // Extra options:
        //
        IntOption    verb   ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 1, IntRange(0, 2));
        IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", INT32_MAX, IntRange(0, INT32_MAX));
        IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", INT32_MAX, IntRange(0, INT32_MAX));
        
        parseOptions(argc, argv, true);

        ERSolver S;
        double initial_time = cpuTime();

        S.verbosity = verb;
        
        solver = &S;
        // Use signal handlers that forcibly quit until the solver will be able to respond to
        // interrupts:
        signal(SIGINT, SIGINT_exit);
        signal(SIGXCPU,SIGINT_exit);

        // Set limit on CPU-time:
        if (cpu_lim != INT32_MAX){
            rlimit rl;
            getrlimit(RLIMIT_CPU, &rl);
            if (rl.rlim_max == RLIM_INFINITY || (rlim_t)cpu_lim < rl.rlim_max){
                rl.rlim_cur = cpu_lim;
                if (setrlimit(RLIMIT_CPU, &rl) == -1)
                    printf("WARNING! Could not set resource limit: CPU-time.\n");
            } }

        // Set limit on virtual memory:
        if (mem_lim != INT32_MAX){
            rlim_t new_mem_lim = (rlim_t)mem_lim * 1024*1024;
            rlimit rl;
            getrlimit(RLIMIT_AS, &rl);
            if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max){
                rl.rlim_cur = new_mem_lim;
                if (setrlimit(RLIMIT_AS, &rl) == -1)
                    printf("WARNING! Could not set resource limit: Virtual memory.\n");
            } }
        
        if (argc == 1)
            printf("Reading from standard input... Use '--help' for help.\n");
        
        gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
        if (in == NULL)
            printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);
        
        if (S.verbosity > 0){
            printf("============================[ Problem Statistics ]=============================\n");
            printf("|                                                                             |\n"); }
        
        parse_DIMACS(in, S);
        gzclose(in);
        FILE* res = (argc >= 3) ? fopen(argv[2], "wb") : NULL;
        
        if (S.verbosity > 0){
            printf("|  Number of variables:  %12d                                         |\n", S.assignmentTrail.nVars());
            printf("|  Number of clauses:    %12d                                         |\n", S.clauseDatabase.nClauses()); }
        
        double parsed_time = cpuTime();
        if (S.verbosity > 0){
            printf("|  Parse time:           %12.2f s                                       |\n", parsed_time - initial_time);
            printf("|                                                                             |\n"); }
 
        // Change to signal-handlers that will only notify the solver and allow it to terminate
        // voluntarily:
        signal(SIGINT, SIGINT_interrupt);
        signal(SIGXCPU,SIGINT_interrupt);
       
        if (!S.simplify()){
            if (res != NULL) fprintf(res, "UNSAT\n"), fclose(res);
            if (S.verbosity > 0){
                printf("===============================================================================\n");
                printf("Solved by unit propagation\n");
                printStats(S);
                printf("\n"); }
            printf("UNSATISFIABLE\n");
            exit(20);
        }
        
        vec<Lit> dummy;
        lbool ret = S.solveLimited(dummy);
        if (S.verbosity > 0){
            printStats(S);
            printf("\n"); }
        printf(ret == l_True ? "SATISFIABLE\n" : ret == l_False ? "UNSATISFIABLE\n" : "INDETERMINATE\n");
        if (res != NULL){
            if (ret == l_True){
                fprintf(res, "SAT\n");
                for (int i = 0; i < S.assignmentTrail.nVars(); i++)
                    if (S.model[i] != l_Undef)
                        fprintf(res, "%s%s%d", (i==0)?"":" ", (S.model[i]==l_True)?"":"-", i+1);
                fprintf(res, " 0\n");
            }else if (ret == l_False)
                fprintf(res, "UNSAT\n");
            else
                fprintf(res, "INDET\n");
            fclose(res);
        }
        
#ifdef NDEBUG
        exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
#else
        return (ret == l_True ? 10 : ret == l_False ? 20 : 0);
#endif
    } catch (OutOfMemoryException&){
        printf("===============================================================================\n");
        printf("INDETERMINATE\n");
        exit(0);
    }
}
