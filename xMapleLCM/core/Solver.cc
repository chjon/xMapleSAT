/***************************************************************************************[Solver.cc]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010, Niklas Sorensson
 
Chanseok Oh's MiniSat Patch Series -- Copyright (c) 2015, Chanseok Oh
 
Maple_LCM, Based on MapleCOMSPS_DRUP -- Copyright (c) 2017, Mao Luo, Chu-Min LI, Fan Xiao: implementing a learnt clause minimisation approach
Reference: M. Luo, C.-M. Li, F. Xiao, F. Manya, and Z. L. , “An effective learnt clause minimization approach for cdcl sat solvers,” in IJCAI-2017, 2017, pp. to–appear.

Maple_LCM_Dist, Based on Maple_LCM -- Copyright (c) 2017, Fan Xiao, Chu-Min LI, Mao Luo: using a new branching heuristic called Distance at the beginning of search

xMaple_LCM_Dist, based on Maple_LCM_Dist -- Copyright (c) 2022, Jonathan Chung, Vijay Ganesh, Sam Buss

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

#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <algorithm>

#include "mtl/Sort.h"
#include "core/Solver.h"
#include "core/SolverER.h"

using namespace Minisat;

#ifdef BIN_DRUP
int Solver::buf_len = 0;
unsigned char Solver::drup_buf[2 * 1024 * 1024];
unsigned char* Solver::buf_ptr = drup_buf;
#endif

//=================================================================================================
// Options:


static const char* _cat = "CORE";

static DoubleOption  opt_clause_decay      (_cat, "cla-decay",   "The clause activity decay factor",              0.999,    DoubleRange(0, false, 1, false));
static IntOption     opt_ccmin_mode        (_cat, "ccmin-mode",  "Controls conflict clause minimization (0=none, 1=basic, 2=deep)", 2, IntRange(0, 2));
static IntOption     opt_restart_first     (_cat, "rfirst",      "The base restart interval", 100, IntRange(1, INT32_MAX));
static DoubleOption  opt_restart_inc       (_cat, "rinc",        "Restart interval increase factor", 2, DoubleRange(1, false, HUGE_VAL, false));
static DoubleOption  opt_garbage_frac      (_cat, "gc-frac",     "The fraction of wasted memory allowed before a garbage collection is triggered",  0.20, DoubleRange(0, false, HUGE_VAL, false));

//=================================================================================================
// Constructor/Destructor:


Solver::Solver() :

    // Parameters (user settable):
    //
    drup_file        (NULL)
  , verbosity        (0)
  , clause_decay     (opt_clause_decay)
  , ccmin_mode       (opt_ccmin_mode)
  , garbage_frac     (opt_garbage_frac)
  , restart_first    (opt_restart_first)
  , restart_inc      (opt_restart_inc)

  // Parameters (the rest):
  //
  , learntsize_factor((double)1/(double)3), learntsize_inc(1.1)

  // Parameters (experimental):
  //
  , learntsize_adjust_start_confl (100)
  , learntsize_adjust_inc         (1.5)


  // Statistics: (formerly in 'SolverStats')
  //
  , solves(0), starts(0), conflicts(0)
  , clauses_literals(0), learnts_literals(0), max_literals(0), tot_literals(0)

  , ok                 (true)
  , cla_inc            (1)
  
  , simpDB_assigns     (-1)
  , simpDB_props       (0)
  , progress_estimate  (0)
  , remove_satisfied   (true)

  , core_lbd_cut       (3)
  , global_lbd_sum     (0)
  , lbd_queue          (50)
  , next_T2_reduce     (10000)
  , next_L_reduce      (15000)
  
  , counter            (0)

  // Resource constraints:
  //
  , conflict_budget    (-1)
  , asynch_interrupt   (false)

  // simplfiy
  , nbSimplifyAll(0)

  // simplifyAll adjust occasion
  , curSimplify(1)
  , nbconfbeforesimplify(1000)
  , incSimplify(1000)

  , branchingHeuristicManager(this)
  , unitPropagator(this)
  , ser (new SolverER(this))
{}

Solver::~Solver() {
    delete ser; ser = nullptr;
}

void Solver::simpleUncheckEnqueue(Lit p, CRef from){
    assert(value(p) == l_Undef);
    variableDatabase.setVar(var(p), lbool(!sign(p)));
    vardata[var(p)].reason = from;
    trail.push_(p);
}

void Solver::cancelUntilTrailRecord() {
    for (int c = trail.size() - 1; c >= trailRecord; c--) {
        Var x = var(trail[c]);
        variableDatabase.setVar(x, l_Undef);
    }

    unitPropagator.setQueueHead(trailRecord);
    trail.shrink(trail.size() - trailRecord);
}

void Solver::litsEnqueue(int cutP, Clause& c) {
    for (int i = cutP; i < c.size(); i++) {
        simpleUncheckEnqueue(~c[i]);
    }
}

bool Solver::removed(CRef cr) {
    return ca[cr].mark() == 1;
}

void Solver::simpleAnalyze(CRef confl, vec<Lit>& out_learnt, vec<CRef>& reason_clause, bool True_confl) {
    int pathC = 0;
    Lit p = lit_Undef;
    int index = trail.size() - 1;

    do{
        if (confl != CRef_Undef){
            reason_clause.push(confl);
            Clause& c = ca[confl];
            // Special case for binary clauses
            // The first one has to be SAT
            if (p != lit_Undef && c.size() == 2 && value(c[0]) == l_False) {
                assert(value(c[1]) == l_True);
                Lit tmp = c[0];
                c[0] = c[1], c[1] = tmp;
            }
            // if True_confl==true, then choose p begin with the 1th index of c;
            for (int j = (p == lit_Undef && True_confl == false) ? 0 : 1; j < c.size(); j++){
                Lit q = c[j];
                if (!seen[var(q)]){
                    seen[var(q)] = 1;
                    pathC++;
                }
            }
        }
        else if (confl == CRef_Undef){
            out_learnt.push(~p);
        }
        // if not break, while() will come to the index of trail blow 0, and fatal error occur;
        if (pathC == 0) break;
        // Select next clause to look at:
        while (!seen[var(trail[index--])]);
        // if the reason cr from the 0-level assigned var, we must break avoid move forth further;
        // but attention that maybe seen[x]=1 and never be clear. However makes no matter;
        if (trailRecord > index + 1) break;
        p = trail[index + 1];
        confl = reason(var(p));
        seen[var(p)] = 0;
        pathC--;

    } while (pathC >= 0);
}

void Solver::simplifyLearnt(Clause& c)
{
    ////
    original_length_record += c.size();

    trailRecord = trail.size();// record the start pointer

    vec<Lit> falseLit;
    falseLit.clear();

    //sort(&c[0], c.size(), VarOrderLevelLt(vardata));

    bool True_confl = false;
    // int beforeSize, afterSize;
    // beforeSize = c.size();
    int i, j;
    CRef confl;

    for (i = 0, j = 0; i < c.size(); i++){
        if (value(c[i]) == l_Undef){
            //printf("///@@@ uncheckedEnqueue:index = %d. l_Undef\n", i);
            simpleUncheckEnqueue(~c[i]);
            c[j++] = c[i];
            confl = unitPropagator.simplePropagate();
            if (confl != CRef_Undef){
                break;
            }
        }
        else{
            if (value(c[i]) == l_True){
                //printf("///@@@ uncheckedEnqueue:index = %d. l_True\n", i);
                c[j++] = c[i];
                True_confl = true;
                confl = reason(var(c[i]));
                break;
            }
            else{
                //printf("///@@@ uncheckedEnqueue:index = %d. l_False\n", i);
                falseLit.push(c[i]);
            }
        }
    }
    c.shrink(c.size() - j);
    // afterSize = c.size();
    //printf("\nbefore : %d, after : %d ", beforeSize, afterSize);


    if (confl != CRef_Undef || True_confl == true){
        simp_learnt_clause.clear();
        simp_reason_clause.clear();
        if (True_confl == true){
            simp_learnt_clause.push(c.last());
        }
        simpleAnalyze(confl, simp_learnt_clause, simp_reason_clause, True_confl);

        if (simp_learnt_clause.size() < c.size()){
            for (i = 0; i < simp_learnt_clause.size(); i++){
                c[i] = simp_learnt_clause[i];
            }
            c.shrink(c.size() - i);
        }
    }

    cancelUntilTrailRecord();

    ////
    simplified_length_record += c.size();

}

bool Solver::simplifyLearnt_x(vec<CRef>& learnts_x)
{
    // int beforeSize, afterSize;
    // int learnts_x_size_before = learnts_x.size();

    int ci, cj, li, lj;
    bool sat, false_lit;
    unsigned int nblevels;
    ////
    //printf("learnts_x size : %d\n", learnts_x.size());

    ////
    int nbSimplified = 0;
    int nbSimplifing = 0;

    for (ci = 0, cj = 0; ci < learnts_x.size(); ci++){
        CRef cr = learnts_x[ci];
        Clause& c = ca[cr];

        if (removed(cr)) continue;
        else if (c.simplified()){
            learnts_x[cj++] = learnts_x[ci];
            ////
            nbSimplified++;
        }
        else{
            ////
            nbSimplifing++;
            sat = false_lit = false;
            for (int i = 0; i < c.size(); i++){
                if (value(c[i]) == l_True){
                    sat = true;
                    break;
                }
                else if (value(c[i]) == l_False){
                    false_lit = true;
                }
            }
            if (sat){
                removeClause(cr);
            }
            else{
                detachClause(cr, true);

                if (false_lit){
                    for (li = lj = 0; li < c.size(); li++){
                        if (value(c[li]) != l_False){
                            c[lj++] = c[li];
                        }
                    }
                    c.shrink(li - lj);
                }

                // beforeSize = c.size();
                assert(c.size() > 1);
                // simplify a learnt clause c
                simplifyLearnt(c);
                assert(c.size() > 0);
                // afterSize = c.size();

                //printf("beforeSize: %2d, afterSize: %2d\n", beforeSize, afterSize);

                if (c.size() == 1){
                    // when unit clause occur, enqueue and propagate
                    uncheckedEnqueue(c[0]);
                    if (unitPropagator.propagate() != CRef_Undef){
                        ok = false;
                        return false;
                    }
                    // delete the clause memory in logic
                    c.mark(1);
                    ca.free(cr);
                }
                else{
                    attachClause(cr);
                    learnts_x[cj++] = learnts_x[ci];

                    nblevels = computeLBD(c);
                    if (nblevels < static_cast<unsigned int>(c.lbd())){
                        //printf("lbd-before: %d, lbd-after: %d\n", c.lbd(), nblevels);
                        c.set_lbd(nblevels);
                    }
                    if (c.mark() != CORE){
                        if (c.lbd() <= core_lbd_cut){
                            //if (c.mark() == LOCAL) local_learnts_dirty = true;
                            //else tier2_learnts_dirty = true;
                            cj--;
                            learnts_core.push(cr);
                            c.mark(CORE);
                        }
                        else if (c.mark() == LOCAL && c.lbd() <= 6){
                            //local_learnts_dirty = true;
                            cj--;
                            learnts_tier2.push(cr);
                            c.mark(TIER2);
                        }
                    }

                    c.setSimplified(true);
                }
            }
        }
    }
    learnts_x.shrink(ci - cj);

    //   printf("c nbLearnts_x %d / %d, nbSimplified: %d, nbSimplifing: %d\n",
    //          learnts_x_size_before, learnts_x.size(), nbSimplified, nbSimplifing);

    return true;
}

bool Solver::simplifyLearnt_core()
{
    // int beforeSize, afterSize;
    // int learnts_core_size_before = learnts_core.size();

    int ci, cj, li, lj;
    bool sat, false_lit;
    unsigned int nblevels;
    ////
    //printf("learnts_x size : %d\n", learnts_x.size());

    ////
    int nbSimplified = 0;
    int nbSimplifing = 0;

    for (ci = 0, cj = 0; ci < learnts_core.size(); ci++){
        CRef cr = learnts_core[ci];
        Clause& c = ca[cr];

        if (removed(cr)) continue;
        else if (c.simplified()){
            learnts_core[cj++] = learnts_core[ci];
            ////
            nbSimplified++;
        }
        else{
            int saved_size=c.size();
            //         if (drup_file){
            //                 add_oc.clear();
            //                 for (int i = 0; i < c.size(); i++) add_oc.push(c[i]); }
            ////
            nbSimplifing++;
            sat = false_lit = false;
            for (int i = 0; i < c.size(); i++){
                if (value(c[i]) == l_True){
                    sat = true;
                    break;
                }
                else if (value(c[i]) == l_False){
                    false_lit = true;
                }
            }
            if (sat){
                removeClause(cr);
            }
            else{
                detachClause(cr, true);

                if (false_lit){
                    for (li = lj = 0; li < c.size(); li++){
                        if (value(c[li]) != l_False){
                            c[lj++] = c[li];
                        }
                    }
                    c.shrink(li - lj);
                }

                // beforeSize = c.size();
                assert(c.size() > 1);
                // simplify a learnt clause c
                simplifyLearnt(c);
                assert(c.size() > 0);
                // afterSize = c.size();
                
                if(drup_file && saved_size !=c.size()){
#ifdef BIN_DRUP
                    binDRUP('a', c , drup_file);
                    //                    binDRUP('d', add_oc, drup_file);
#else
                    for (int i = 0; i < c.size(); i++)
                        fprintf(drup_file, "%i ", (var(c[i]) + 1) * (-2 * sign(c[i]) + 1));
                    fprintf(drup_file, "0\n");

                    //                    fprintf(drup_file, "d ");
                    //                    for (int i = 0; i < add_oc.size(); i++)
                    //                        fprintf(drup_file, "%i ", (var(add_oc[i]) + 1) * (-2 * sign(add_oc[i]) + 1));
                    //                    fprintf(drup_file, "0\n");
#endif
                }

                //printf("beforeSize: %2d, afterSize: %2d\n", beforeSize, afterSize);

                if (c.size() == 1){
                    // when unit clause occur, enqueue and propagate
                    uncheckedEnqueue(c[0]);
                    if (unitPropagator.propagate() != CRef_Undef){
                        ok = false;
                        return false;
                    }
                    // delete the clause memory in logic
                    c.mark(1);
                    ca.free(cr);
//#ifdef BIN_DRUP
//                    binDRUP('d', c, drup_file);
//#else
//                    fprintf(drup_file, "d ");
//                    for (int i = 0; i < c.size(); i++)
//                        fprintf(drup_file, "%i ", (var(c[i]) + 1) * (-2 * sign(c[i]) + 1));
//                    fprintf(drup_file, "0\n");
//#endif
                }
                else{
                    attachClause(cr);
                    learnts_core[cj++] = learnts_core[ci];

                    nblevels = computeLBD(c);
                    if (nblevels < static_cast<unsigned int>(c.lbd())){
                        //printf("lbd-before: %d, lbd-after: %d\n", c.lbd(), nblevels);
                        c.set_lbd(nblevels);
                    }

                    c.setSimplified(true);
                }
            }
        }
    }
    learnts_core.shrink(ci - cj);

    //    printf("c nbLearnts_core %d / %d, nbSimplified: %d, nbSimplifing: %d\n",
    //           learnts_core_size_before, learnts_core.size(), nbSimplified, nbSimplifing);

    return true;

}

bool Solver::simplifyLearnt_tier2()
{
    // int beforeSize, afterSize;
    // int learnts_tier2_size_before = learnts_tier2.size();

    int ci, cj, li, lj;
    bool sat, false_lit;
    unsigned int nblevels;
    ////
    //printf("learnts_x size : %d\n", learnts_x.size());

    ////
    int nbSimplified = 0;
    int nbSimplifing = 0;

    for (ci = 0, cj = 0; ci < learnts_tier2.size(); ci++){
        CRef cr = learnts_tier2[ci];
        Clause& c = ca[cr];

        if (removed(cr)) continue;
        else if (c.simplified()){
            learnts_tier2[cj++] = learnts_tier2[ci];
            ////
            nbSimplified++;
        }
        else{
            int saved_size=c.size();
            //            if (drup_file){
            //                    add_oc.clear();
            //                    for (int i = 0; i < c.size(); i++) add_oc.push(c[i]); }
            ////
            nbSimplifing++;
            sat = false_lit = false;
            for (int i = 0; i < c.size(); i++){
                if (value(c[i]) == l_True){
                    sat = true;
                    break;
                }
                else if (value(c[i]) == l_False){
                    false_lit = true;
                }
            }
            if (sat){
                removeClause(cr);
            }
            else{
                detachClause(cr, true);

                if (false_lit){
                    for (li = lj = 0; li < c.size(); li++){
                        if (value(c[li]) != l_False){
                            c[lj++] = c[li];
                        }
                    }
                    c.shrink(li - lj);
                }

                // beforeSize = c.size();
                assert(c.size() > 1);
                // simplify a learnt clause c
                simplifyLearnt(c);
                assert(c.size() > 0);
                // afterSize = c.size();
                
                if(drup_file && saved_size!=c.size()){

#ifdef BIN_DRUP
                    binDRUP('a', c , drup_file);
                    //                    binDRUP('d', add_oc, drup_file);
#else
                    for (int i = 0; i < c.size(); i++)
                        fprintf(drup_file, "%i ", (var(c[i]) + 1) * (-2 * sign(c[i]) + 1));
                    fprintf(drup_file, "0\n");

                    //                    fprintf(drup_file, "d ");
                    //                    for (int i = 0; i < add_oc.size(); i++)
                    //                        fprintf(drup_file, "%i ", (var(add_oc[i]) + 1) * (-2 * sign(add_oc[i]) + 1));
                    //                    fprintf(drup_file, "0\n");
#endif
                }

                //printf("beforeSize: %2d, afterSize: %2d\n", beforeSize, afterSize);

                if (c.size() == 1){
                    // when unit clause occur, enqueue and propagate
                    uncheckedEnqueue(c[0]);
                    if (unitPropagator.propagate() != CRef_Undef){
                        ok = false;
                        return false;
                    }
                    // delete the clause memory in logic
                    c.mark(1);
                    ca.free(cr);
//#ifdef BIN_DRUP
//                    binDRUP('d', c, drup_file);
//#else
//                    fprintf(drup_file, "d ");
//                    for (int i = 0; i < c.size(); i++)
//                        fprintf(drup_file, "%i ", (var(c[i]) + 1) * (-2 * sign(c[i]) + 1));
//                    fprintf(drup_file, "0\n");
//#endif
                }
                else{
                    attachClause(cr);
                    learnts_tier2[cj++] = learnts_tier2[ci];

                    nblevels = computeLBD(c);
                    if (nblevels < static_cast<unsigned int>(c.lbd())){
                        //printf("lbd-before: %d, lbd-after: %d\n", c.lbd(), nblevels);
                        c.set_lbd(nblevels);
                    }

                    if (c.lbd() <= core_lbd_cut){
                        cj--;
                        learnts_core.push(cr);
                        c.mark(CORE);
                    }

                    c.setSimplified(true);
                }
            }
        }
    }
    learnts_tier2.shrink(ci - cj);

    //    printf("c nbLearnts_tier2 %d / %d, nbSimplified: %d, nbSimplifing: %d\n",
    //           learnts_tier2_size_before, learnts_tier2.size(), nbSimplified, nbSimplifing);

    return true;

}

bool Solver::simplifyAll()
{
    ////
    simplified_length_record = original_length_record = 0;

    if (!ok || unitPropagator.propagate() != CRef_Undef)
        return ok = false;

    //// cleanLearnts(also can delete these code), here just for analyzing
    //if (local_learnts_dirty) cleanLearnts(learnts_local, LOCAL);
    //if (tier2_learnts_dirty) cleanLearnts(learnts_tier2, TIER2);
    //local_learnts_dirty = tier2_learnts_dirty = false;

    if (!simplifyLearnt_core()) return ok = false;
    if (!simplifyLearnt_tier2()) return ok = false;
    //if (!simplifyLearnt_x(learnts_local)) return ok = false;

    checkGarbage();

    ////
    //  printf("c size_reduce_ratio     : %4.2f%%\n",
    //         original_length_record == 0 ? 0 : (original_length_record - simplified_length_record) * 100 / (double)original_length_record);

    return true;
}
//=================================================================================================
// Minor methods:


// Creates a new SAT variable in the solver. If 'decision' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool sign, bool dvar) {
    int v = variableDatabase.newVar();
    vardata.push(mkVarData(CRef_Undef, 0));

    seen .push(0);
    seen2.push(0);
    trail.capacity(v+1);

    branchingHeuristicManager  .newVar(v, sign, dvar);
    unitPropagator.newVar(v);
    ser->newVar();

    return v;
}


bool Solver::addClause_(vec<Lit>& ps)
{
    assert(decisionLevel() == 0);
    if (!ok) return false;

    // Check if clause is satisfied and remove false/duplicate literals:
    sort(ps);
    Lit p; int i, j;

    if (drup_file){
        add_oc.clear();
        for (int i = 0; i < ps.size(); i++) add_oc.push(ps[i]); }

    for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
        if (value(ps[i]) == l_True || ps[i] == ~p)
            return true;
        else if (value(ps[i]) != l_False && ps[i] != p)
            ps[j++] = p = ps[i];
    ps.shrink(i - j);

    if (drup_file && i != j){
#ifdef BIN_DRUP
        binDRUP('a', ps, drup_file);
        binDRUP('d', add_oc, drup_file);
#else
        for (int i = 0; i < ps.size(); i++)
            fprintf(drup_file, "%i ", (var(ps[i]) + 1) * (-2 * sign(ps[i]) + 1));
        fprintf(drup_file, "0\n");

        fprintf(drup_file, "d ");
        for (int i = 0; i < add_oc.size(); i++)
            fprintf(drup_file, "%i ", (var(add_oc[i]) + 1) * (-2 * sign(add_oc[i]) + 1));
        fprintf(drup_file, "0\n");
#endif
    }

    if (ps.size() == 0)
        return ok = false;
    else if (ps.size() == 1){
        uncheckedEnqueue(ps[0]);
        return ok = (unitPropagator.propagate() == CRef_Undef);
    }else{
        CRef cr = ca.alloc(ps, false);
        clauses.push(cr);
        attachClause(cr);
    }

    return true;
}


void Solver::attachClause(CRef cr) {
    const Clause& c = ca[cr];
    unitPropagator.attachClause(c, cr);
    if (c.learnt()) learnts_literals += c.size();
    else            clauses_literals += c.size();
}


void Solver::detachClause(CRef cr, bool strict) {
    const Clause& c = ca[cr];
    unitPropagator.detachClause(c, cr, strict);
    if (c.learnt()) learnts_literals -= c.size();
    else            clauses_literals -= c.size();
}


void Solver::removeClause(CRef cr) {
    Clause& c = ca[cr];

    if (drup_file){
        if (c.mark() != 1){
#ifdef BIN_DRUP
            binDRUP('d', c, drup_file);
#else
            fprintf(drup_file, "d ");
            for (int i = 0; i < c.size(); i++)
                fprintf(drup_file, "%i ", (var(c[i]) + 1) * (-2 * sign(c[i]) + 1));
            fprintf(drup_file, "0\n");
#endif
        }else
            printf("c Bug. I don't expect this to happen.\n");
    }

    detachClause(cr);
    // Don't leave pointers to free'd memory!
    if (locked(c)){
        Lit implied = c.size() != 2 ? c[0] : (value(c[0]) == l_True ? c[0] : c[1]);
        vardata[var(implied)].reason = CRef_Undef; }
    c.mark(1);
    ca.free(cr);
}


bool Solver::satisfied(const Clause& c) const {
    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) == l_True)
            return true;
    return false; }


// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void Solver::cancelUntil(int level) {
    if (decisionLevel() <= level) return;

    // Unassign variables set later than the given level
    for (int c = trail.size()-1; c >= trail_lim[level]; c--){
        Var      x  = var(trail[c]);
        variableDatabase.setVar(x, l_Undef);
        branchingHeuristicManager.handleEventLitUnassigned(trail[c], conflicts, c > trail_lim.last());
    }

    // Update the propagation queue
    unitPropagator.setQueueHead(trail_lim[level]);

    // Update the assignment trail
    trail.shrink(trail.size() - trail_lim[level]);
    trail_lim.shrink(trail_lim.size() - level);
}


//=================================================================================================
// Major methods:


/*_________________________________________________________________________________________________
|
|  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
|  
|  Description:
|    Analyze conflict and produce a reason clause.
|  
|    Pre-conditions:
|      * 'out_learnt' is assumed to be cleared.
|      * Current decision level must be greater than root level.
|  
|    Post-conditions:
|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
|      * If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of the 
|        rest of literals. There may be others from the same level though.
|  
|________________________________________________________________________________________________@*/
void Solver::analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel, int& out_lbd)
{
    int pathC = 0;
    Lit p     = lit_Undef;

    // Generate conflict clause:
    //
    out_learnt.push();      // (leave room for the asserting literal)
    int index   = trail.size() - 1;

    do {
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

        // For binary clauses, we don't rearrange literals in propagate(), so check and make sure the first is an implied lit.
        if (p != lit_Undef && c.size() == 2 && value(c[0]) == l_False){
            assert(value(c[1]) == l_True);
            std::swap(c[0], c[1]);
        }

        // Update LBD if improved.
        if (c.learnt() && c.mark() != CORE){
            int lbd = computeLBD(c);
            if (lbd < c.lbd()){
                if (c.lbd() <= 30) c.removable(false); // Protect once from reduction.
                c.set_lbd(lbd);
                if (lbd <= core_lbd_cut){
                    learnts_core.push(confl);
                    c.mark(CORE);
                }else if (lbd <= 6 && c.mark() == LOCAL){
                    // Bug: 'cr' may already be in 'learnts_tier2', e.g., if 'cr' was demoted from TIER2
                    // to LOCAL previously and if that 'cr' is not cleaned from 'learnts_tier2' yet.
                    learnts_tier2.push(confl);
                    c.mark(TIER2); }
            }

            if (c.mark() == TIER2)
                c.touched() = conflicts;
            else if (c.mark() == LOCAL)
                claBumpActivity(c);
        }

        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];

            if (!seen[var(q)] && level(var(q)) > 0) {
                branchingHeuristicManager.handleEventLitInConflictGraph(q);
                seen[var(q)] = 1;
                if (level(var(q)) >= decisionLevel()) {
                    pathC++;
                } else {
                    out_learnt.push(q);
                }
            }
        }
        
        // Select next clause to look at:
        while (!seen[var(trail[index--])]);
        p     = trail[index+1];
        confl = reason(var(p));
        seen[var(p)] = 0;
        pathC--;

    } while (pathC > 0);
    out_learnt[0] = ~p;

    // Simplify conflict clause:
    //
    int i, j;
    out_learnt.copyTo(analyze_toclear);
    if (ccmin_mode == 2){
        uint32_t abstract_level = 0;
        for (i = 1; i < out_learnt.size(); i++)
            abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)

        for (i = j = 1; i < out_learnt.size(); i++)
            if (reason(var(out_learnt[i])) == CRef_Undef || !litRedundant(out_learnt[i], abstract_level))
                out_learnt[j++] = out_learnt[i];
        
    } else if (ccmin_mode == 1) {
        for (i = j = 1; i < out_learnt.size(); i++){
            Var x = var(out_learnt[i]);

            if (reason(x) == CRef_Undef)
                out_learnt[j++] = out_learnt[i];
            else{
                Clause& c = ca[reason(var(out_learnt[i]))];
                for (int k = c.size() == 2 ? 0 : 1; k < c.size(); k++)
                    if (!seen[var(c[k])] && level(var(c[k])) > 0){
                        out_learnt[j++] = out_learnt[i];
                        break; }
            }
        }
    }else {
        i = j = out_learnt.size();
    }

    max_literals += out_learnt.size();
    out_learnt.shrink(i - j);
    tot_literals += out_learnt.size();

    out_lbd = computeLBD(out_learnt);
    if (out_lbd <= 6 && out_learnt.size() <= 30) // Try further minimization?
        if (binResMinimize(out_learnt))
            out_lbd = computeLBD(out_learnt); // Recompute LBD if minimized.

    // Find correct backtrack level:
    //
    if (out_learnt.size() == 1) {
        out_btlevel = 0;
    } else {
        int max_i = 1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 2; i < out_learnt.size(); i++)
            if (level(var(out_learnt[i])) > level(var(out_learnt[max_i])))
                max_i = i;
        // Swap-in this literal at index 1:
        Lit p             = out_learnt[max_i];
        out_learnt[max_i] = out_learnt[1];
        out_learnt[1]     = p;
        out_btlevel       = level(var(p));
    }
}


// Try further learnt clause minimization by means of binary clause resolution.
bool Solver::binResMinimize(vec<Lit>& out_learnt) {
    // Preparation: remember which false variables we have in 'out_learnt'.
    counter++;
    for (int i = 1; i < out_learnt.size(); i++)
        seen2[var(out_learnt[i])] = counter;

    // Get the list of binary clauses containing 'out_learnt[0]'.
    const vec<Watcher>& ws = unitPropagator.getBinaryWatchers(~out_learnt[0]);

    int to_remove = 0;
    for (int i = 0; i < ws.size(); i++){
        Lit the_other = ws[i].blocker;
        // Does 'the_other' appear negatively in 'out_learnt'?
        if (seen2[var(the_other)] == counter && value(the_other) == l_True){
            to_remove++;
            seen2[var(the_other)] = counter - 1; // Remember to remove this variable.
        }
    }

    // Shrink.
    if (to_remove > 0){
        int last = out_learnt.size() - 1;
        for (int i = 1; i < out_learnt.size() - to_remove; i++)
            if (seen2[var(out_learnt[i])] != counter)
                out_learnt[i--] = out_learnt[last--];
        out_learnt.shrink(to_remove);
    }
    return to_remove != 0;
}


// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, uint32_t abstract_levels) {
    analyze_stack.clear(); analyze_stack.push(p);
    int top = analyze_toclear.size();
    while (analyze_stack.size() > 0){
        assert(reason(var(analyze_stack.last())) != CRef_Undef);
        Clause& c = ca[reason(var(analyze_stack.last()))]; analyze_stack.pop();

        // Special handling for binary clauses like in 'analyze()'.
        if (c.size() == 2 && value(c[0]) == l_False){
            assert(value(c[1]) == l_True);
            Lit tmp = c[0];
            c[0] = c[1], c[1] = tmp; }

        for (int i = 1; i < c.size(); i++){
            Lit p  = c[i];
            if (!seen[var(p)] && level(var(p)) > 0){
                if (reason(var(p)) != CRef_Undef && (abstractLevel(var(p)) & abstract_levels) != 0){
                    seen[var(p)] = 1;
                    analyze_stack.push(p);
                    analyze_toclear.push(p);
                }else{
                    for (int j = top; j < analyze_toclear.size(); j++)
                        seen[var(analyze_toclear[j])] = 0;
                    analyze_toclear.shrink(analyze_toclear.size() - top);
                    return false;
                }
            }
        }
    }

    return true;
}


/*_________________________________________________________________________________________________
|
|  analyzeFinal : (p : Lit)  ->  [void]
|  
|  Description:
|    Specialized analysis procedure to express the final conflict in terms of assumptions.
|    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
|    stores the result in 'out_conflict'.
|________________________________________________________________________________________________@*/
void Solver::analyzeFinal(Lit p, vec<Lit>& out_conflict) {
    out_conflict.clear();
    out_conflict.push(p);

    if (decisionLevel() == 0)
        return;

    seen[var(p)] = 1;

    for (int i = trail.size()-1; i >= trail_lim[0]; i--){
        Var x = var(trail[i]);
        if (seen[x]){
            if (reason(x) == CRef_Undef){
                assert(level(x) > 0);
                out_conflict.push(~trail[i]);
            }else{
                Clause& c = ca[reason(x)];
                for (int j = c.size() == 2 ? 0 : 1; j < c.size(); j++)
                    if (level(var(c[j])) > 0)
                        seen[var(c[j])] = 1;
            }
            seen[x] = 0;
        }
    }

    seen[var(p)] = 0;
}


void Solver::uncheckedEnqueue(Lit p, CRef from) {
    assert(value(p) == l_Undef);
    Var x = var(p);
    variableDatabase.setVar(x, lbool(!sign(p)));
    branchingHeuristicManager.handleEventLitAssigned(p, conflicts); // Update heuristic
    vardata[x] = mkVarData(from, decisionLevel());
    trail.push_(p);
}

/*_________________________________________________________________________________________________
|
|  reduceDB : ()  ->  [void]
|  
|  Description:
|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
|________________________________________________________________________________________________@*/
struct reduceDB_lt { 
    ClauseAllocator& ca;
    reduceDB_lt(ClauseAllocator& ca_) : ca(ca_) {}
    bool operator () (CRef x, CRef y) const { return ca[x].activity() < ca[y].activity(); }
};
void Solver::reduceDB() {
    int     i, j;
    //if (local_learnts_dirty) cleanLearnts(learnts_local, LOCAL);
    //local_learnts_dirty = false;

    sort(learnts_local, reduceDB_lt(ca));

    int limit = learnts_local.size() / 2;
    for (i = j = 0; i < learnts_local.size(); i++){
        Clause& c = ca[learnts_local[i]];
        if (c.mark() == LOCAL)
            if (c.removable() && !locked(c) && i < limit)
                removeClause(learnts_local[i]);
            else{
                if (!c.removable()) limit++;
                c.removable(true);
                learnts_local[j++] = learnts_local[i]; }
    }
    learnts_local.shrink(i - j);

    checkGarbage();
}
void Solver::reduceDB_Tier2() {
    int i, j;
    for (i = j = 0; i < learnts_tier2.size(); i++){
        Clause& c = ca[learnts_tier2[i]];
        if (c.mark() == TIER2)
            if (!locked(c) && c.touched() + 30000 < conflicts){
                learnts_local.push(learnts_tier2[i]);
                c.mark(LOCAL);
                //c.removable(true);
                c.activity() = 0;
                claBumpActivity(c);
            }else
                learnts_tier2[j++] = learnts_tier2[i];
    }
    learnts_tier2.shrink(i - j);
}


void Solver::removeSatisfied(vec<CRef>& cs) {
    int i, j;
    for (i = j = 0; i < cs.size(); i++){
        Clause& c = ca[cs[i]];
        if (satisfied(c))
            removeClause(cs[i]);
        else
            cs[j++] = cs[i];
    }
    cs.shrink(i - j);
}

void Solver::safeRemoveSatisfied(vec<CRef>& cs, unsigned valid_mark) {
    int i, j;
    for (i = j = 0; i < cs.size(); i++){
        Clause& c = ca[cs[i]];
        if (c.mark() == valid_mark)
            if (satisfied(c))
                removeClause(cs[i]);
            else
                cs[j++] = cs[i];
    }
    cs.shrink(i - j);
}


/*_________________________________________________________________________________________________
|
|  simplify : [void]  ->  [bool]
|  
|  Description:
|    Simplify the clause database according to the current top-level assigment. Currently, the only
|    thing done here is the removal of satisfied clauses, but more things can be put here.
|________________________________________________________________________________________________@*/
bool Solver::simplify()
{
    assert(decisionLevel() == 0);

    if (!ok || unitPropagator.propagate() != CRef_Undef)
        return ok = false;

    if (nAssigns() == simpDB_assigns || (simpDB_props > 0))
        return true;

    // Remove satisfied clauses:
    removeSatisfied(learnts_core); // Should clean core first.
    safeRemoveSatisfied(learnts_tier2, TIER2);
    safeRemoveSatisfied(learnts_local, LOCAL);
    if (remove_satisfied) {        // Can be turned off.
        removeSatisfied(clauses);
        ser->removeSatisfied();
    }
    checkGarbage();
    branchingHeuristicManager.rebuildOrderHeap();

    simpDB_assigns = nAssigns();
    simpDB_props   = clauses_literals + learnts_literals;   // (shouldn't depend on stats really, but it will do for now)

    return true;
}

CRef Solver::propagateLits(vec<Lit>& lits) {
    Lit lit;
    int i;

    for(i=lits.size()-1; i>=0; i--) {
        lit=lits[i];
        if (value(lit) == l_Undef) {
            newDecisionLevel();
            uncheckedEnqueue(lit);
            CRef confl = unitPropagator.propagate();
            if (confl != CRef_Undef) {
                return confl;
            }
        }
    }
    return CRef_Undef;
}
/*_________________________________________________________________________________________________
|
|  search : (nof_conflicts : int) (params : const SearchParams&)  ->  [lbool]
|  
|  Description:
|    Search for a model the specified number of conflicts. 
|  
|  Output:
|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
|________________________________________________________________________________________________@*/
lbool Solver::search(int& nof_conflicts) {
    assert(ok);
    int         backtrack_level;
    int         lbd;
    vec<Lit>    learnt_clause;
    bool        cached = false;
    starts++;

    // simplify
    //
    if (conflicts >= static_cast<uint64_t>(curSimplify * nbconfbeforesimplify)){
        //        printf("c ### simplifyAll on conflict : %lld\n", conflicts);
        //printf("nbClauses: %d, nbLearnts_core: %d, nbLearnts_tier2: %d, nbLearnts_local: %d, nbLearnts: %d\n",
        //	clauses.size(), learnts_core.size(), learnts_tier2.size(), learnts_local.size(),
        //	learnts_core.size() + learnts_tier2.size() + learnts_local.size());
        nbSimplifyAll++;
        if (!simplifyAll()){
            return l_False;
        }
        curSimplify = (conflicts / nbconfbeforesimplify) + 1;
        nbconfbeforesimplify += incSimplify;
    }

#if ER_USER_GEN_LOCATION == ER_GEN_LOCATION_AFTER_RESTART
    // Generate extension variable definitions
    // Only try generating more extension variables if there aren't any buffered already
    ser->generateDefinitions();
#endif

#if ER_USER_ADD_LOCATION == ER_ADD_LOCATION_AFTER_RESTART
    // Add extension variables
    ser->introduceExtVars(ser->extDefs);
#endif

    for (;;){
        CRef confl = unitPropagator.propagate();

        if (confl != CRef_Undef){
            // CONFLICT
            conflicts++; nof_conflicts--;
            branchingHeuristicManager.handleEventConflicted(conflicts);

            if (conflicts == 100000 && learnts_core.size() < 100) core_lbd_cut = 5;
            if (decisionLevel() == 0) return l_False;

            if (branchingHeuristicManager.currentBranchingHeuristic() == BranchingHeuristic::VSIDS &&
                branchingHeuristicManager.distanceBranchingEnabled()
            ) {
                branchingHeuristicManager.collectFirstUIP(confl);
            }

            // Generate learnt clause
            learnt_clause.clear();
            analyze(confl, learnt_clause, backtrack_level, lbd);
            branchingHeuristicManager.handleEventLearnedClause(learnt_clause, backtrack_level);

            cancelUntil(backtrack_level);

            // EXTENDED RESOLUTION - substitute disjunctions with extension variables
            // This must be called after backtracking because extension variables might need to be propagated
#if ER_USER_SUBSTITUTE_HEURISTIC != ER_SUBSTITUTE_HEURISTIC_NONE
            ser->substitute(learnt_clause, ser->user_extSubPredicate);
#endif

            lbd--;
            if (branchingHeuristicManager.currentBranchingHeuristic() == BranchingHeuristic::VSIDS){
                cached = false;
                lbd_queue.push(lbd);
                global_lbd_sum += (lbd > 50 ? 50 : lbd);
            }

            // Add learnt clause to clause database
            if (learnt_clause.size() == 1){
                uncheckedEnqueue(learnt_clause[0]);
            }else{
                CRef cr = ca.alloc(learnt_clause, true);
                ca[cr].set_lbd(lbd);
                if (lbd <= core_lbd_cut){
                    learnts_core.push(cr);
                    ca[cr].mark(CORE);
                }else if (lbd <= 6){
                    learnts_tier2.push(cr);
                    ca[cr].mark(TIER2);
                    ca[cr].touched() = conflicts;
                }else{
                    learnts_local.push(cr);
                    claBumpActivity(ca[cr]); }
                attachClause(cr);
                uncheckedEnqueue(learnt_clause[0], cr);

#if ER_ENABLE_GLUCOSER
                static FilterPredicate fil_ler = std::bind(&SolverER::user_extFilPredicate_ler, ser, std::placeholders::_1);
                ser->filterIncremental(cr, fil_ler, SolverER::HeuristicType::LER);
#endif
#if ER_USER_FILTER_HEURISTIC != ER_FILTER_HEURISTIC_NONE
                ser->filterIncremental(cr, ser->user_extFilPredicate, SolverER::HeuristicType::DEFAULT);
#endif
            }

            // Output learnt clause to proof file
            if (drup_file){
#ifdef BIN_DRUP
                binDRUP('a', learnt_clause, drup_file);
#else
                for (int i = 0; i < learnt_clause.size(); i++)
                    fprintf(drup_file, "%i ", (var(learnt_clause[i]) + 1) * (-2 * sign(learnt_clause[i]) + 1));
                fprintf(drup_file, "0\n");
#endif
            }

            // Update activities
            branchingHeuristicManager.decayActivityVSIDS();
            claDecayActivity();

            /*if (--learntsize_adjust_cnt == 0){
                learntsize_adjust_confl *= learntsize_adjust_inc;
                learntsize_adjust_cnt    = (int)learntsize_adjust_confl;
                max_learnts             *= learntsize_inc;

                if (verbosity >= 1)
                    printf("c | %9d | %7d %8d %8d | %8d %8d %6.0f | %6.3f %% |\n",
                           (int)conflicts,
                           (int)branchingHeuristicManager.dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(), (int)clauses_literals,
                           (int)max_learnts, nLearnts(), (double)learnts_literals/nLearnts(), progressEstimate()*100);
            }*/

#if ER_ENABLE_GLUCOSER
            {
                using namespace std::placeholders;

                // Generate extension variable definitions
                // Only try generating more extension variables if there aren't any buffered already
                static SelectionHeuristic sel_ler = std::bind(&SolverER::user_extSelHeuristic_all, ser, _1, _2, _3);
                static ExtDefHeuristic    def_ler = std::bind(&SolverER::user_extDefHeuristic_ler, ser, _1, _2, _3);
                ser->selectClauses(sel_ler, SolverER::HeuristicType::LER, 1);
                ser->defineExtVars(def_ler, SolverER::HeuristicType::LER);

                // Add extension variables
                ser->introduceExtVars(ser->extDefs, SolverER::HeuristicType::LER);
            }
#endif

        } else {
            // NO CONFLICT
            bool restart = false;
            if (branchingHeuristicManager.currentBranchingHeuristic() == BranchingHeuristic::CHB)
                restart = nof_conflicts <= 0;
            else if (!cached){
                restart = lbd_queue.full() && (lbd_queue.avg() * 0.8 > global_lbd_sum / branchingHeuristicManager.conflicts_VSIDS);
                cached = true;
            }
            if (restart || !withinBudget()){
                lbd_queue.clear();
                cached = false;
                // Reached bound on number of conflicts:
                progress_estimate = progressEstimate();
                cancelUntil(0);
                return l_Undef; }

            // Simplify the set of problem clauses:
            if (decisionLevel() == 0) {
                if (!simplify()) return l_False;
#if ER_USER_DELETE_HEURISTIC != ER_DELETE_HEURISTIC_NONE
                ser->deleteExtVarsIfNecessary();
#endif
            }

            if (conflicts >= next_T2_reduce){
                next_T2_reduce = conflicts + 10000;
                reduceDB_Tier2(); }
            if (conflicts >= next_L_reduce){
                next_L_reduce = conflicts + 15000;
                reduceDB(); }

            Lit next = lit_Undef;
            {
                // New variable decision:
                next = branchingHeuristicManager.pickBranchLit();
                if (next == lit_Undef)
                    // Model found:
                    return l_True;

                if (ser->isExtVar(var(next))) ser->branchOnExt++;
            }

            // Increase decision level and enqueue 'next'
            newDecisionLevel();
            uncheckedEnqueue(next);
        }
    }
}


double Solver::progressEstimate() const {
    double  progress = 0;
    double  F = 1.0 / variableDatabase.nVars();

    for (int i = 0; i <= decisionLevel(); i++){
        int beg = i == 0 ? 0 : trail_lim[i - 1];
        int end = i == decisionLevel() ? trail.size() : trail_lim[i];
        progress += pow(F, i) * (end - beg);
    }

    return progress / variableDatabase.nVars();
}

/*
  Finite subsequences of the Luby-sequence:

  0: 1
  1: 1 1 2
  2: 1 1 2 1 1 2 4
  3: 1 1 2 1 1 2 4 1 1 2 1 1 2 4 8
  ...


 */

static double luby(double y, int x){

    // Find the finite subsequence that contains index 'x', and the
    // size of that subsequence:
    int size, seq;
    for (size = 1, seq = 0; size < x+1; seq++, size = 2*size+1);

    while (size-1 != x){
        size = (size-1)>>1;
        seq--;
        x = x % size;
    }

    return pow(y, seq);
}

// NOTE: assumptions passed in member-variable 'assumptions'.
lbool Solver::solve_()
{
    //signal(SIGALRM, SIGALRM_switch);
    //alarm(2500);

    model.clear();
    conflict.clear();
    if (!ok) return l_False;

    solves++;

    max_learnts               = nClauses() * learntsize_factor;
    learntsize_adjust_confl   = learntsize_adjust_start_confl;
    learntsize_adjust_cnt     = (int)learntsize_adjust_confl;
    lbool   status            = l_Undef;

    if (verbosity >= 1){
        printf("c ============================[ Search Statistics ]==============================\n");
        printf("c | Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
        printf("c |           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
        printf("c ===============================================================================\n");
    }

    add_tmp.clear();

    ser->originalNumVars = variableDatabase.nVars();

    // Warmup period: search with VSIDS
    branchingHeuristicManager.setBranchingHeuristic(BranchingHeuristic::VSIDS);
    int init = 10000;
    while (status == l_Undef && init > 0 && withinBudget())
        status = search(init);
    branchingHeuristicManager.setBranchingHeuristic(BranchingHeuristic::CHB);

    // Search:
    int curr_restarts = 0;
    while (status == l_Undef && withinBudget()) {
        // Periodically switch branching heuristic
        branchingHeuristicManager.handleEventRestarted(unitPropagator.propagations);
        if (branchingHeuristicManager.currentBranchingHeuristic() == BranchingHeuristic::VSIDS) {
            int weighted = INT32_MAX;
            status = search(weighted);
        } else {
            int nof_conflicts = luby(restart_inc, curr_restarts) * restart_first;
            curr_restarts++;
            status = search(nof_conflicts);
        }
    }

    if (verbosity >= 1)
        printf("c ===============================================================================\n");

#ifdef BIN_DRUP
    if (drup_file && status == l_False) binDRUP_flush(drup_file);
#endif

    if (status == l_True){
        // Extend & copy model:
        model.growTo(variableDatabase.nVars());
        for (int i = 0; i < variableDatabase.nVars(); i++) model[i] = value(i);
    }else if (status == l_False && conflict.size() == 0)
        ok = false;

    cancelUntil(0);
    return status;
}

//=================================================================================================
// Writing CNF to DIMACS:
// 
// FIXME: this needs to be rewritten completely.

static Var mapVar(Var x, vec<Var>& map, Var& max)
{
    if (map.size() <= x || map[x] == -1){
        map.growTo(x+1, -1);
        map[x] = max++;
    }
    return map[x];
}


void Solver::toDimacs(FILE* f, Clause& c, vec<Var>& map, Var& max)
{
    if (satisfied(c)) return;

    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) != l_False)
            fprintf(f, "%s%d ", sign(c[i]) ? "-" : "", mapVar(var(c[i]), map, max)+1);
    fprintf(f, "0\n");
}


void Solver::toDimacs(const char *file, const vec<Lit>& assumps)
{
    FILE* f = fopen(file, "wr");
    if (f == NULL)
        fprintf(stderr, "could not open file %s\n", file), exit(1);
    toDimacs(f, assumps);
    fclose(f);
}


void Solver::toDimacs(FILE* f, const vec<Lit>& assumps)
{
    // Handle case when solver is in contradictory state:
    if (!ok){
        fprintf(f, "p cnf 1 2\n1 0\n-1 0\n");
        return; }

    vec<Var> map; Var max = 0;

    // Cannot use removeClauses here because it is not safe
    // to deallocate them at this point. Could be improved.
    int cnt = 0;
    for (int i = 0; i < clauses.size(); i++)
        if (!satisfied(ca[clauses[i]]))
            cnt++;

    for (int i = 0; i < clauses.size(); i++)
        if (!satisfied(ca[clauses[i]])){
            Clause& c = ca[clauses[i]];
            for (int j = 0; j < c.size(); j++)
                if (value(c[j]) != l_False)
                    mapVar(var(c[j]), map, max);
        }

    // Assumptions are added as unit clauses:
    cnt += assumptions.size();

    fprintf(f, "p cnf %d %d\n", max, cnt);

    for (int i = 0; i < assumptions.size(); i++){
        assert(value(assumptions[i]) != l_False);
        fprintf(f, "%s%d 0\n", sign(assumptions[i]) ? "-" : "", mapVar(var(assumptions[i]), map, max)+1);
    }

    for (int i = 0; i < clauses.size(); i++)
        toDimacs(f, ca[clauses[i]], map, max);

    if (verbosity > 0)
        printf("c Wrote %d clauses with %d variables.\n", cnt, max);
}


//=================================================================================================
// Garbage Collection methods:

void Solver::relocAll(ClauseAllocator& to)
{
    // All watchers:
    //
    unitPropagator.relocAll(to);

    // All reasons:
    //
    for (int i = 0; i < trail.size(); i++){
        Var v = var(trail[i]);

        if (reason(v) != CRef_Undef && (ca[reason(v)].reloced() || locked(ca[reason(v)])))
            ca.reloc(vardata[v].reason, to);
    }

    // All learnt:
    //
    for (int i = 0; i < learnts_core.size(); i++)
        ca.reloc(learnts_core[i], to);
    for (int i = 0; i < learnts_tier2.size(); i++)
        ca.reloc(learnts_tier2[i], to);
    for (int i = 0; i < learnts_local.size(); i++)
        ca.reloc(learnts_local[i], to);

    // All extension variable definition clauses
    //
    ser->relocAll(to);

    // All original:
    //
    int i, j;
    for (i = j = 0; i < clauses.size(); i++)
        if (ca[clauses[i]].mark() != 1){
            ca.reloc(clauses[i], to);
            clauses[j++] = clauses[i]; }
    clauses.shrink(i - j);
}


void Solver::garbageCollect()
{
    // Initialize the next region to a size corresponding to the estimated utilization degree. This
    // is not precise but should avoid some unnecessary reallocations for the new region:
    ClauseAllocator to(ca.size() - ca.wasted());

    relocAll(to);
    if (verbosity >= 2)
        printf("c |  Garbage collection:   %12d bytes => %12d bytes             |\n",
               ca.size()*ClauseAllocator::Unit_Size, to.size()*ClauseAllocator::Unit_Size);
    to.moveTo(ca);
}
