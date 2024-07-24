/****************************************************************************************[Solver.C]
 GlucoeER 1.0 -- Copyright (c) 2009-2010, Gilles Audemard, George Katsirelos, Laurent Simon
                                CRIL - Univ. Artois, France
                                LRI  - Univ. Paris Sud, France

 Glucose 1.0 -- Copyright (c) 2009, Gilles Audemard, Laurent Simon
                                CRIL - Univ. Artois, France
                                LRI  - Univ. Paris Sud, France

GlucosER is based on MiniSat. The licence of Glucose is exactly the
same as Minisat (see below)

 MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson

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


#include "Solver.h"
#include "Sort.h"
#include <cmath>
#include "Constants.h"
#include <iostream>
#include <iomanip>
using namespace std;

#ifdef PRINTDEBUG
#define DOUT if(0); else cout
#else
#define DOUT if(0) cout
#endif

#ifdef PRINTDEBUG
#define DBG if(0); else
#else
#define DBG if(0)
#endif


struct p_lit {
  Lit l;
  p_lit(Lit ll) : l(ll) {}
};

Clause* Clause_new(const vec<Lit>& ps, bool learnt = false) {
  assert(sizeof(Lit)      == sizeof(uint32_t));
  assert(sizeof(float)    == sizeof(uint32_t));
  void* mem = malloc(sizeof(Clause) + sizeof(uint32_t)*(ps.size()));
  return new (mem) Clause(ps, learnt); }


std::ostream& operator<<(std::ostream& os, p_lit const& p)
{
  os << (sign(p.l) ? "-" : "") << (var(p.l)+1);
  return os;
}

template<typename clausetype>
struct clause_printer {
  clausetype& c;
  clause_printer(clausetype& cc) : c(cc) {}
};

template<typename T>
std::ostream& operator<<(std::ostream& os, clause_printer<T> const& p)
{
  os << '(';
  for(int i = 0; i != p.c.size(); ++i)
    os << setw(8) << p_lit(p.c[i]) << ((i+1) % 8 == 0 ? '\n' : ' ');
  os << ')';
  return os;
}

template<typename clausetype>
struct clause_lvl_printer {
  clausetype& c;
  Solver *s;
  clause_lvl_printer(clausetype& cc, Solver* ss) : c(cc), s(ss) {}
};

template<typename T>
std::ostream& operator<<(std::ostream& os, clause_lvl_printer<T> const& p)
{
  os << "<\n";
  for(int i = 0; i != p.c.size(); ++i)
    os << setw(8) << p_lit(p.c[i])
       << " ("
       << ((p.s->value(p.c[i]) == l_Undef) ? -1 : p.s->level[var(p.c[i])])
       << ")"
       << ((i+1) % 8 == 0 ? "\n" : " ");
  os << "\n>\n";
  return os;
}

clause_printer<Clause> print_clause(Clause* c)
{
  return clause_printer<Clause>(*c);
}

clause_printer<vec<Lit> > print_veclit(vec<Lit>& c)
{
  return clause_printer<vec<Lit> >(c);
}

clause_lvl_printer<Clause> print_clause_lvl(Clause* c, Solver *s)
{
  return clause_lvl_printer<Clause>(*c, s);
}

clause_lvl_printer<vec<Lit> > print_clause_lvl(vec<Lit>& c, Solver *s)
{
  return clause_lvl_printer<vec<Lit> >(c, s);
}

std::ostream& operator<<(std::ostream& os, lbool l)
{
  if( l == l_True )
    os << "TRUE";
  else if( l == l_False )
    os << "FALSE";
  else if( l == l_Undef )
    os << "undef";
  else
    assert(0);
  return os;
}


//================================================================================================
// Function to detach binary and other clauses

static inline void removeW(vec<Watched> &ws,Clause *c) {
  int j = 0;
  for (; j < ws.size() && ws[j].wcl != c; j++);
  assert(j < ws.size());
  for (; j < ws.size()-1; j++) ws[j] = ws[j+1];
  ws.pop();

}

static inline void removeBin(vec<Binaire> &ws,Clause *c) {
  int j = 0;
  for (; j < ws.size() && ws[j].clause != c; j++);
  assert(j < ws.size());
  for (; j < ws.size()-1; j++) ws[j] = ws[j+1];
  ws.pop();

}


//=================================================================================================
// Constructor/Destructor:





Solver::Solver() :

    // Parameters: (formerly in 'SearchParams')
    var_decay(1 / 0.95), clause_decay(1 / 0.999), random_var_freq(0.02)
  , restart_first(100), restart_inc(1.5), learntsize_factor((double)1/(double)3), learntsize_inc(1)

    // More parameters:
    //
  , expensive_ccmin  (true)
  , polarity_mode    (polarity_user)
  , verbosity        (0)

    // Statistics: (formerly in 'SolverStats')
    //
    , nbDL2(0),nbBin(0),nbUn(0) , nbReduceDB(0), starts(0), decisions(0), rnd_decisions(0), propagations(0), conflicts(0)
  , clauses_literals(0), learnts_literals(0), max_literals(0), tot_literals(0)
  , extensions(0), extensions_noninput(0)
  , lit_replacements(0)
  , extend_branch(0)
  , ext_neg_lits(0), ext_pos_lits(0)

  , ok               (true)
  , cla_inc          (1)
  , var_inc          (1)
  , qhead            (0)
  , simpDB_assigns   (-1)
  , simpDB_props     (0)
  , order_heap       (VarOrderLt(activity))
  , random_seed      (91648253)
  , progress_estimate(0)
  , remove_satisfied (true)
    , nbclausesbeforereduce (NBCLAUSESBEFOREREDUCE)
    , extend_reduce_candidates( compare_vsids(activity) )

{MYFLAG = 0;}


Solver::~Solver()
{
    for (int i = 0; i < learnts.size(); i++) free(learnts[i]);
    for (int i = 0; i < clauses.size(); i++) free(clauses[i]);
}


//=================================================================================================
// Minor methods:


// Creates a new SAT variable in the solver. If 'decision_var' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool sign, bool dvar)
{
    int v = nVars();
    watches   .push();          // (list for positive literal)
    watches   .push();          // (list for negative literal)
    watchesBin   .push();          // (binary clauses list for positive literal)
    watchesBin   .push();          // (binary clauses list for negative literal)
    reason    .push(NULL);
    assigns   .push(toInt(l_Undef));
    level     .push(-1);
    activity  .push(0);
    seen      .push(0);
    permDiff  .push(0);
    permDiff  .push(0);

    extended  .push(); // 1 entry for each literal
    extended  .push();
    extending .push();
    extending_clauses.push();
    extended_blocked_count.push();
    extended_blocks.push();
    extend_deleted.push();
    extend_reduce_index.push();


    polarity    .push((char)sign);
    decision_var.push((char)dvar);

    insertVarOrder(v);
    return v;
}


bool Solver::addClause(vec<Lit>& ps)
{
  assert(decisionLevel() == 0);

  if (!ok)
    return false;
  else{
    // Check if clause is satisfied and remove false/duplicate literals:
    sort(ps);
        Lit p; int i, j;
        for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
          if (value(ps[i]) == l_True || ps[i] == ~p)
            return true;
          else if (value(ps[i]) != l_False && ps[i] != p)
            ps[j++] = p = ps[i];
        ps.shrink(i - j);
  }

  if (ps.size() == 0)
    return ok = false;
  else if (ps.size() == 1){
    assert(value(ps[0]) == l_Undef);
    uncheckedEnqueue(ps[0]);
    return ok = (propagate() == NULL);
  }else{
    Clause* c = Clause_new(ps, false);
    clauses.push(c);
    attachClause(*c);
  }

  return true;
}


void Solver::attachClause(Clause& c) {
  assert(c.size() > 1);
  if(c.size()==2) {
    watchesBin[toInt(~c[0])].push();
    watchesBin[toInt(~c[1])].push();
    watchesBin[toInt(~c[0])].last().clause = &c;
    watchesBin[toInt(~c[0])].last().implied = c[1];
    watchesBin[toInt(~c[1])].last().clause = &c;
    watchesBin[toInt(~c[1])].last().implied = c[0];

  } else {
    watches[toInt(~c[0])].push();
    watches[toInt(~c[1])].push();
    watches[toInt(~c[0])].last().wcl = &c;
    watches[toInt(~c[0])].last().blocked = c[c.size()/2];
    watches[toInt(~c[1])].last().wcl = &c;
    watches[toInt(~c[1])].last().blocked = c[c.size()/2];
  }
  if (c.learnt()) learnts_literals += c.size();
  else            clauses_literals += c.size();
}

void Solver::detachClause(Clause& c) {
  assert(c.size() > 1);
  //    assert(find(watches[toInt(~c[0])], &c));
  //assert(find(watches[toInt(~c[1])], &c));

  if(c.size()==2) {
    removeBin(watchesBin[toInt(~c[0])],&c);
    removeBin(watchesBin[toInt(~c[1])],&c);
  } else {
    removeW(watches[toInt(~c[0])], &c);
    removeW(watches[toInt(~c[1])], &c);
  }
  if (c.learnt()) learnts_literals -= c.size();
  else            clauses_literals -= c.size(); }


void Solver::removeClause(Clause& c) {
  detachClause(c);
  free(&c); }


bool Solver::satisfied(const Clause& c) const {
  for (int i = 0; i < c.size(); i++)
    if (value(c[i]) == l_True || extend_deleted[var(c[i])])
      return true;
  return false; }

// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void Solver::cancelUntil(int level) {
  if (decisionLevel() > level){
    for (int c = trail.size()-1; c >= trail_lim[level]; c--){
      Var     x  = var(trail[c]);
      assigns[x] = toInt(l_Undef);
      insertVarOrder(x); }
    qhead = trail_lim[level];
    trail.shrink(trail.size() - trail_lim[level]);
    trail_lim.shrink(trail_lim.size() - level);
  }
}


//=================================================================================================
// Major methods:


Lit Solver::pickBranchLit(int polarity_mode, double random_var_freq)
{
  Var next = var_Undef;

  // Random decision:
  if ((drand(random_seed) < random_var_freq) && !order_heap.empty()){

    next = order_heap[irand(random_seed,order_heap.size())];
    if (toLbool(assigns[next]) == l_Undef && decision_var[next])
      rnd_decisions++; }

  // Activity based decision:
  while (next == var_Undef || toLbool(assigns[next]) != l_Undef || !decision_var[next]
         || extend_deleted[next])
    if (order_heap.empty()){
      next = var_Undef;
      break;
    }else
      next = order_heap.removeMin();
  if( is_extended(next) ) ++ extend_branch;

  bool sign = false;
  switch (polarity_mode){
  case polarity_true:  sign = false; break;
  case polarity_false: sign = true;  break;
  case polarity_user:  if(next!=var_Undef) sign = polarity[next]; break;
  case polarity_rnd:   sign = irand(random_seed, 2); break;
  default: assert(false); }

  return next == var_Undef ? lit_Undef : Lit(next, sign);
}


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
|
|  Effect:
|    Will undo part of the trail, upto but not beyond the assumption of the current decision level.
|________________________________________________________________________________________________@*/
void Solver::analyze(Clause* confl, vec<Lit>& out_learnt, int& out_btlevel,int &lbd)
{
  int pathC = 0;
  Lit p     = lit_Undef;


  out_learnt.push();      // (leave room for the asserting literal)
  int index   = trail.size() - 1;
  out_btlevel = 0;

  do{
    assert(confl != NULL);          // (otherwise should be UIP)
    Clause& c = *confl;

    // Special case for binary clauses
    // The first one has to be SAT
    if( p != lit_Undef && c.size()==2 && value(c[0])==l_False) {
      assert(value(c[1])==l_True);
      Lit tmp = c[0];
      c[0] =  c[1], c[1] = tmp;
    }

    if (c.learnt())
      claBumpActivity(c);

    for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
      Lit q = c[j];

      if (!seen[var(q)] && level[var(q)] > 0){

        varBumpActivity(var(q));

        seen[var(q)] = 1;
        if (level[var(q)] >= decisionLevel()){
          pathC++;
#ifdef UPDATEVARACTIVITY
          // UPDATEVARACTIVITY trick (see competition'09 companion paper)
          if((reason[var(q)]!=NULL)  && (reason[var(q)]->learnt()))
            lastDecisionLevel.push(q);
#endif
        }
        else{
          out_learnt.push(q);
          if (level[var(q)] > out_btlevel)
            out_btlevel = level[var(q)];
        }
      }
    }
    // Select next clause to look at:
    while (!seen[var(trail[index--])]);
    p     = trail[index+1];
    confl = reason[var(p)];
    seen[var(p)] = 0;
    pathC--;

  }while (pathC > 0);
  out_learnt[0] = ~p;



    // Simplify conflict clause:
    //
  int i, j;
  out_learnt.copyTo(analyze_toclear);



  if (expensive_ccmin){
    uint32_t abstract_level = 0;
    for (i = 1; i < out_learnt.size(); i++)
      abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)

    for (i = j = 1; i < out_learnt.size(); i++)
      if (reason[var(out_learnt[i])] == NULL || !litRedundant(out_learnt[i], abstract_level))
        out_learnt[j++] = out_learnt[i];
  }else{
    /*    for (i = j = 1; i < out_learnt.size(); i++){
      Clause& c = *reason[var(out_learnt[i])];
      if(c.size()==2 && value(c[0])==l_False) {
        assert(value(c[1])==l_True);
        Lit tmp = c[0];
        c[0] =  c[1], c[1] = tmp;
      }

      for (int k = 1; k < c.size(); k++)
        if (!seen[var(c[k])] && level[var(c[k])] > 0){
          out_learnt[j++] = out_learnt[i];
          break; }
          }*/
    i = j = 0;
  }
  max_literals += out_learnt.size();
  out_learnt.shrink(i - j);
  tot_literals += out_learnt.size();




  extend_reduce_assertLits.clear();
  extend_reduce_assertClauses.clear();

  j = extend_reduce(out_learnt);
  out_learnt.shrink(j);


  // ALBERT
  // for(int i = 0; i != out_learnt.size(); ++i) {
  //   if( is_extended(var(out_learnt[i])) ) {
  //     if( sign(out_learnt[i]) )
  //       ++ext_neg_lits;
  //     else
  //       ++ext_pos_lits;
  //   }
  // }


    int extend_btlevel;

    // Not here
    //j = extend_new_reduce(out_learnt, extend_btlevel);
    j = 0;
    //j = extend_level_reduce(out_learnt, extend_btlevel);
    out_learnt.shrink(j);
    //DOUT << "\tnew extensions shrank it by another " << j << " to \n"
    //     << print_veclit(out_learnt) << "\n";

    bool assert_head = false;

    tot_literals += out_learnt.size();

    // Find correct backtrack level:
    //
    if (out_learnt.size() == 1)
        out_btlevel = 0;
    else{
      //if( assert_head ) {
      int max_i = 1;
      int second_watch = -1;
      while( value(out_learnt[max_i]) == l_Undef &&
             max_i < out_learnt.size() ) {
        assert_head = false;
        second_watch = max_i;
        ++max_i;
      }
      for (int i = max_i+1; i < out_learnt.size(); i++) {
        if( value(out_learnt[i]) == l_Undef ) {
          assert_head = false;
          second_watch = i;
        }
        if (level[var(out_learnt[i])] > level[var(out_learnt[max_i])])
          max_i = i;
      }
      Lit p             = out_learnt[max_i];
      out_btlevel       = level[var(p)];
      if( second_watch > 0 ) {
        swap(out_learnt[1], out_learnt[second_watch]);
      } else {
        out_learnt[max_i] = out_learnt[1];
        out_learnt[1]     = p;
      }

      if( second_watch > 0 )
        out_btlevel = extend_btlevel;
    }






  // Find the LBD measure
  lbd = 0;
  MYFLAG++;
  for(int i=0;i<out_learnt.size();i++) {

    int l = level[var(out_learnt[i])];
    if (permDiff[l] != MYFLAG) {
      permDiff[l] = MYFLAG;
      lbd++;
    }
  }



#ifdef UPDATEVARACTIVITY
  // UPDATEVARACTIVITY trick (see competition'09 companion paper)
  if(lastDecisionLevel.size()>0) {
    for(int i = 0;i<lastDecisionLevel.size();i++) {
      if(reason[var(lastDecisionLevel[i])]->lbd()<lbd)
        varBumpActivity(var(lastDecisionLevel[i]));
    }
    lastDecisionLevel.clear();
  }
#endif




    for (int j = 0; j < analyze_toclear.size(); j++) seen[var(analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)

    //for (int i = 0; i != nVars(); ++i ) assert(seen[i] == 0);

}



/*____________________________________________________________
| extend_reduce : (out_learnt : vec<Lit>) -> int
|
| Find pairs of literals that have already been used to extend
| the formula and use them to shrink this clause
|_____________________________________________________________*/

// this uses bit 32 of seen[] to signal that a variable has been
// replaced by an extended variable and bit 64 to signal that a
// variable is in the clause
int Solver::extend_reduce(vec<Lit>& out_learnt)
{
  // ALBERT
  return 0;
  static const int replbit = 32;
  static const int inclbit = 64;
  //  DOUT << "Reducing clause " << print_clause_lvl(out_learnt, this) << "\n";

  int osize = out_learnt.size();
  Heap< compare_vsids >& candidates = extend_reduce_candidates;
  assert( candidates.empty() );

  // we start at 1: whatever we do, we leave the asserting literal alone
  for(int i = 1; i != osize; ++i) {
    extend_reduce_index[var(out_learnt[i])] = i;
    seen[var(out_learnt[i])] |= inclbit;
    for(int j = i+1; j != osize; ++j) {
      Lit l1 = out_learnt[i], l2 = out_learnt[j];
      if( toInt(l1) > toInt(l2) )
        swap(l1, l2);
      //      //DOUT << "Checking whether X <=> " << p_lit(l1) << "\\/" << p_lit(l2)
      //<< " has been introduced\n";
      extmap& el1 = extended[toInt(l1)];
      extmap::iterator emi = el1.find(toInt(l2));
      if( emi == el1.end() ) continue;
      int X = emi->second;
      if(X) {
        candidates.insert( X );
      }
    }
  }

  while( !candidates.empty() ) {
    int X = candidates.removeMin();
    Lit lX(X);
    Lit l1 = extending[X].first;
    Lit l2 = extending[X].second;
    if( (seen[var(l1)] & replbit) || (seen[var(l2)] & replbit) )
      continue;

    assert(value(lX) == l_False);
    // if( value(lX) != l_False ) {
    //   extend_reduce_assertLits.push(lX);
    //   extend_reduce_assertClauses.push(MPIRI);
    // }

    int idx1 = extend_reduce_index[var(l1)],
      idx2 = extend_reduce_index[var(l2)];
    assert( out_learnt[idx1] == l1 );
    assert( out_learnt[idx2] == l2 );
    if( idx1 > idx2 )
      swap(idx1, idx2);

    assert(idx1 < osize);
    assert(idx2 < osize);

    if( seen[X] & inclbit ) {
      // X is already in the clause, just remove l1 and l2

      out_learnt[ idx2 ] = out_learnt[osize-2];
      out_learnt[ idx1 ] = out_learnt[osize-1];
      extend_reduce_index[var(out_learnt[idx2])] = idx2;
      extend_reduce_index[var(out_learnt[idx1])] = idx1;
      osize -= 2;
    } else {
      out_learnt[ idx1 ] = lX;
      out_learnt[ idx2 ] = out_learnt[osize-1];
      extend_reduce_index[X] = idx1;
      extend_reduce_index[var(out_learnt[osize-1])] = idx2;
      --osize;
    }
    if( !seen[var(l1)] )
      analyze_toclear.push(l1);
    if( !seen[var(l2)] )
      analyze_toclear.push(l2);
    seen[var(l1)] |= replbit;
    seen[var(l2)] |= replbit;
    if( !seen[X] ) {
      seen[X] = 1 | inclbit;
      analyze_toclear.push(lX);
    }

    ++lit_replacements;


    // now add to the heap new candidates that extend Y <=> X \/ l3
    if( !seen[X] ) {
      int i = 1;
      for(int j = i+1; j != osize; ++j) {
        Lit l1 = lX, l2 = out_learnt[j];
        if( toInt(l1) > toInt(l2) )
          swap(l1, l2);
        extmap& el1 = extended[toInt(l1)];
        extmap::iterator emi = el1.find(toInt(l2));
        if( emi == el1.end() ) continue;
        int X = emi->second;
        if(X)
          candidates.insert( X );
      }
    }
  }

  return out_learnt.size() - osize;
}



/* Introduce X <=> lx1 \/ lx2 */
int Solver::extend_by(Lit lx1, Lit lx2, Clause **assertClause)
{
  // ALBERT
  return -1;
  if( toInt(lx2) < toInt(lx1) )
    swap(lx1, lx2);

  {
    extmap& extl1 = extended[toInt(lx1)];
    if( extl1.find(toInt(lx2)) != extl1.end() )
      return -1;
  }

  ++extensions;
  if( is_extended(var(lx1)) || is_extended(var(lx2)) )
    ++extensions_noninput;

  int X;
  if( extend_freelist.size() ) {
    X = extend_freelist.last();
    extend_freelist.pop();
    extend_deleted[X] = false;
  } else
    X = newVar(true,true);
  Lit lX(X);

  extendvars.push(X);
  if( is_extended(var(lx1)) ) {
    ++extended_blocked_count[var(lx1)];
    extended_blocks[X].push(var(lx1));
  }
  if( is_extended(var(lx2)) ) {
    ++extended_blocked_count[var(lx2)];
    extended_blocks[X].push(var(lx2));
  }


  extmap& extl1 = extended[toInt(lx1)];

  extl1[toInt(lx2)] = X;
  extending[X] = make_pair(lx1, lx2);

  vec<Lit> v1, v2, v3;
  // X <-> x1 \/ x2
  v1.push(~lX);
  if( value(lx1) == l_Undef ||
      (value(lx2) == l_False &&
       level[var(lx1)] >= level[var(lx2)] )) {
    v1.push(lx1);
    v1.push(lx2);
  } else {
    v1.push(lx2);
    v1.push(lx1);
  }

  v2.push(lX);
  v2.push(~lx1);


  v3.push(lX);
  v3.push(~lx2);


  Clause* c1 = Clause_new(v1, false);
  //clauses.push(c1);
  extending_clauses[X].xl1l2 = c1;
  attachClause(*c1);
  *assertClause = c1;
  if(value(lx1)==l_False && value(lx2)==l_False) {
    uncheckedEnqueue(~lX,c1);
    //      assert((*assertClause)[1]==~l2);
  }


  Clause* c3 = Clause_new(v3, false);
  //clauses.push(c3);
  extending_clauses[X].xl2 = c3;
  attachClause(*c3);

  if(value(lx2)==l_True) {
    uncheckedEnqueue(lX,c3);
  }


  Clause* c2 = Clause_new(v2, false);
  //clauses.push(c2);
  extending_clauses[X].xl1 = c2;
  attachClause(*c2);

 if(value(lx1)==l_True && value(lx2)!=l_True) {
    uncheckedEnqueue(lX,c2);
  }


  return X;
}


void Solver::deleteVar(int X) {
  // ALBERt
  return;
  assert(value(X) == l_Undef);

  pair<Lit, Lit> ex = extending[X];
  extending[X] = make_pair(lit_Undef, lit_Undef);
  extend_deleted[X] = true;
  extend_freelist.push(X);

  Lit lx1 = ex.first, lx2 = ex.second;
  extmap& extl1 = extended[toInt(lx1)];
  extl1.erase(toInt(lx2));

  // remove x <=> lx1 v lx2 clauses
  extclauses& xc = extending_clauses[X];
  removeClause( *xc.xl1l2 );
  removeClause( *xc.xl1 );
  removeClause( *xc.xl2 );
  xc.xl1l2 = xc.xl1 = xc.xl2 = NULL;

  for(int i = 0; i != extended_blocks[X].size(); ++i) {
    int Xp = extended_blocks[X][i];
    --extended_blocked_count[Xp];
  }
  extended_blocks[X].clear();
}

bool Solver::tryToExtend2(vec<Lit>& cur,int curDBL) {
  return false;
}


/**********************************************
  tryToExtend

  last is the precedent Assertive clause
  cur  is the current   Assertive clause

  We look for
  last = l1 v A v B
  cur  = l2 v A v C

  We do not want to have |l1| in cur !

  intersect will become A

*/


bool Solver::tryToExtend(vec<Lit>& last,int lastDBL,vec<Lit>& cur,int curDBL,vec<Lit>& intersect) {

  bool taut = false;
  if(last.size()==0) return false; // At the begining of the search
  int nb = 0;
  Lit l1 = last[0];
  Lit l2 = cur[0];
  if(var(l1)==var(l2)) return false;
  //printf("   -- try : ");printLit(l1);printf("  ");printLit(l2);printf("\n");
  assert(value(l2)==l_True);
  //assert(value(l1)!=l_False);
  //cout << "--\n" << print_clause_lvl(cur, this) << "\n";
  //cout << "--\n" << print_clause_lvl(last, this) << "\n";

  MYFLAG++; // init

  // Current Clause
  for(int i = 1;i<cur.size();i++) {
    permDiff[toInt(cur[i])] = MYFLAG;
  }


  // |l1| not in cur
  if(permDiff[toInt(l1)]==MYFLAG ||  permDiff[toInt(~l1)]==MYFLAG){
    return false;
  }

  // Last Clause
  for(int i = 1;i<last.size();i++) {
    if(permDiff[toInt(last[i])] == MYFLAG) {
      nb++;//intersect.push(last[i]);
    }

    if(permDiff[toInt(~last[i])]== MYFLAG) {
      printf("TAUT \n");
      ///printLit(last[i]);
      return false; // A litteral and opposite are in the clause
    }
  }


  if(permDiff[toInt(l2)]==MYFLAG ||  permDiff[toInt(~l2)]==MYFLAG){
    printf("HOP\n");
    return false;
  }

  if(nb<=2) return false;
  //cout << "--\n" << print_clause_lvl(cur, this) << "\n";
  //cout << "--\n" << print_clause_lvl(last, this) << "\n";
  //exit(1);

  //  if(!(last.size()-1-nb==0 || cur.size()-1-nb==0)) return false;

  if(last.size()-1-nb!=0 || cur.size()-1-nb!=0) return false; // both C and B are empty


  //  if((lastDBL+curDBL>6)) return false;

  //  if(lastDBL>2 && curDBL>2) return false;


  // Now we can extend the formula
  //  printf("YES\n");
  // cout << "--\n" << print_clause_lvl(intersect, this) << "\n";

  //  cout << "--\n" << print_clause_lvl(cur, this) << "\n";
  //cout << "--\n" << print_clause_lvl(last, this) << "\n";

  Clause *assertClause = 0L;


  int X = extend_by(~l1,~l2, &assertClause);

    for(int i = 0;i<lastExtended.size();i++) {
    extend_by(~l2,lastExtended[i],&assertClause);
  }

  lastExtended.push(~l1);

  if(X==-1) return false;
  Lit lX(X);
  return true;
  if(1 || last.size()-1-nb==0 || cur.size()-1-nb==0) {
    //if(last.size()+cur.size()-2-nb<10) {// ATTENTINO
    // Create the new clause ~X v A v B v C
    vec<Lit> added;
    added.push(~lX);
    for(int i = 1;i<cur.size();i++) {
      added.push(cur[i]);
    }
    for(int i = 1;i<last.size();i++) {
      if(permDiff[toInt(last[i])] != MYFLAG)
        added.push(last[i]);
    }
    Clause* c1 = Clause_new(added, true);
    //clauses.push(c1);
    learnts.push(c1);
    attachClause(*c1);
    c1->setLBD(lastDBL+curDBL-2);
  }



    return true;



}


// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, uint32_t abstract_levels)
{
  analyze_stack.clear(); analyze_stack.push(p);
  int top = analyze_toclear.size();
  while (analyze_stack.size() > 0){
    assert(reason[var(analyze_stack.last())] != NULL);
    Clause& c = *reason[var(analyze_stack.last())]; analyze_stack.pop();
    if(c.size()==2 && value(c[0])==l_False) {
      assert(value(c[1])==l_True);
      Lit tmp = c[0];
      c[0] =  c[1], c[1] = tmp;
    }


    for (int i = 1; i < c.size(); i++){
      Lit p  = c[i];
      if (!seen[var(p)] && level[var(p)] > 0){
        if (reason[var(p)] != NULL && (abstractLevel(var(p)) & abstract_levels) != 0){
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
void Solver::analyzeFinal(Lit p, vec<Lit>& out_conflict)
{
  out_conflict.clear();
  out_conflict.push(p);

  if (decisionLevel() == 0)
    return;

  seen[var(p)] = 1;

  for (int i = trail.size()-1; i >= trail_lim[0]; i--){
    Var x = var(trail[i]);
    if (seen[x]){
      if (reason[x] == NULL){
        assert(level[x] > 0);
        out_conflict.push(~trail[i]);
      }else{
        Clause& c = *reason[x];
        for (int j = 1; j < c.size(); j++)
          if (level[var(c[j])] > 0)
            seen[var(c[j])] = 1;
      }
      seen[x] = 0;
    }
  }

  seen[var(p)] = 0;
}


void Solver::uncheckedEnqueue(Lit p, Clause* from)
{

  assert(value(p) == l_Undef);
  assigns [var(p)] = toInt(lbool(!sign(p)));  // <<== abstract but not uttermost effecient
  level   [var(p)] = decisionLevel();
  reason  [var(p)] = from;
  polarity[var(p)] = sign(p); // We always use phase caching  !
  trail.push(p);
}


/*_________________________________________________________________________________________________
|
|  propagate : [void]  ->  [Clause*]
|
|  Description:
|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
|    otherwise NULL.
|
|    Post-conditions:
|      * the propagation queue is empty, even if there was a conflict.
|________________________________________________________________________________________________@*/
Clause* Solver::propagate()
{
  Clause* confl     = NULL;
  int     num_props = 0;

  while (qhead < trail.size()){
    Lit            p   = trail[qhead++];     // 'p' is enqueued fact to propagate.
    vec<Watched>&  ws  = watches[toInt(p)];
    Watched         *i, *j, *end;
    num_props++;


    // First, Propagate binary clauses
    vec<Binaire> & wbin = watchesBin[toInt(p)];

    for(int k = 0;k<wbin.size();k++) {
      Lit imp = wbin[k].implied;
      if(value(imp) == l_False) {
        return wbin[k].clause;
      }

      if(value(imp) == l_Undef) {
        assert( !extend_deleted[var(imp)] );
        uncheckedEnqueue(imp,wbin[k].clause);
      }
    }


    for (i = j = (Watched*)ws, end = i + ws.size();  i != end;){
      if(value(i->blocked)==l_True) { // Clause is sat
        *j++ = *i++;
        continue;
      }
      Lit bl = i->blocked;
      Clause& c = *(i->wcl);
      i++;


      // Make sure the false literal is data[1]:
      Lit false_lit = ~p;
      if (c[0] == false_lit)
        c[0] = c[1], c[1] = false_lit;

      assert(c[1] == false_lit);

      // If 0th watch is true, then clause is already satisfied.
      Lit first = c[0];
      if (value(first) == l_True){
        j->wcl = &c;
        j->blocked = first;
        j++;
      }else{
        // Look for new watch:
        for (int k = 2; k < c.size(); k++)
          if (value(c[k]) != l_False){
            c[1] = c[k]; c[k] = false_lit;
            watches[toInt(~c[1])].push();
            watches[toInt(~c[1])].last().wcl = &c;
            watches[toInt(~c[1])].last().blocked = c[0];
            goto FoundWatch; }

        // Did not find watch -- clause is unit under assignment:
        j->wcl = &c;
        j->blocked = bl;
        j++;

        if (value(first) == l_False){
          confl = &c;
          qhead = trail.size();
          // Copy the remaining watches:
          while (i < end)
            *j++ = *i++;
        }else {
          assert( !extend_deleted[var(first)] );
          uncheckedEnqueue(first, &c);

#ifdef DYNAMICNBLEVEL
          // DYNAMIC NBLEVEL trick (see competition'09 companion paper)
          if(c.learnt()  && c.lbd()>2) {
            MYFLAG++;
            int nblevels =0;
            for(int i=0;i<c.size();i++) {
              int l = level[var(c[i])];
              if (permDiff[l] != MYFLAG) {
                permDiff[l] = MYFLAG;
                nblevels++;
              }


            }
            if(nblevels+1<c.lbd()) {
              c.setLBD(nblevels);
            }
          }
#endif
        }
      }
    FoundWatch:;
    }
    ws.shrink(i - j);
  }
  propagations += num_props;
  simpDB_props -= num_props;

  return confl;
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
  bool operator () (Clause* x, Clause* y) {

    // Main criteria... Like in MiniSat we keep all binary clauses
    if(x->size()> 2 && y->size()==2) return 1;

    if(y->size()>2 && x->size()==2) return 0;
    if(x->size()==2 && y->size()==2) return 0;

    // Second one  based on literal block distance
    if(x->lbd()> y->lbd()) return 1;
    if(x->lbd()< y->lbd()) return 0;


    // Finally we can use old activity or size, we choose the last one
    //    return x->oldActivity() < y->oldActivity();
    return x->size() < y->size();
  }};




void Solver::reduceDB()
{
  int     i, j;


  nbReduceDB++;
  sort(learnts, reduceDB_lt());

  for (i = j = 0; i < learnts.size() / RATIOREMOVECLAUSES; i++){
    if (learnts[i]->size() > 2 && !locked(*learnts[i]) && learnts[i]->lbd()>2){
      removeClause(*learnts[i]);
    }
    else
      learnts[j++] = learnts[i];
  }
  for (; i < learnts.size(); i++){
    learnts[j++] = learnts[i];
  }
  learnts.shrink(i - j);

  schedule_reduce_extend = true;
}


void Solver::removeSatisfied(vec<Clause*>& cs)
{
  int i,j;
  for (i = j = 0; i < cs.size(); i++){
    if (satisfied(*cs[i]))
      removeClause(*cs[i]);
    else
      cs[j++] = cs[i];
  }
  cs.shrink(i - j);
}

/*_________________________________________________________________________________________________
|
|  reduceExtended : ()  ->  [void]
|
|  Description:
|    Remove extended variables.
|________________________________________________________________________________________________@*/
void Solver::reduceExtended()
{
  if( !schedule_reduce_extend )
    return;

  sort(extendvars, InvVarOrderLt(activity));

  int i, j;
  for (i = j = 0; i < extendvars.size() / RATIOREMOVEVARS; i++){
    int X = extendvars[i];
    if( extended_blocked_count[X] == 0 &&
        value(X) == l_Undef )
      deleteVar(X);
    else
      extendvars[j++] = extendvars[i];
  }
  for (; i < extendvars.size(); i++){
    extendvars[j++] = extendvars[i];
  }
  extendvars.shrink(i - j);
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

  if (!ok || propagate() != NULL)
    return ok = false;

  if ( !schedule_reduce_extend &&
       (nAssigns() == simpDB_assigns || (simpDB_props > 0))) {
    return true;
  }

  schedule_reduce_extend = false;

  // Remove satisfied clauses:
  removeSatisfied(learnts);
  if (remove_satisfied)        // Can be turned off.
    removeSatisfied(clauses);

  // Remove fixed variables from the variable heap:
  order_heap.filter(VarFilter(*this));

  simpDB_assigns = nAssigns();
  simpDB_props   = clauses_literals + learnts_literals;   // (shouldn't depend on stats really, but it will do for now)

  return true;
}


/*_________________________________________________________________________________________________
|
|  search : ()
|
|  Description:
|    Search for a model the specified number of conflicts, keeping the number of learnt clauses
|    below the provided limit.
|
|  Output:
|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.

|________________________________________________________________________________________________@*/

#define LAST 3

static long curRestart = 1,cons=0,conf4Stats=0;

lbool Solver::search()
{
  assert(ok);

  int         backtrack_level;
  int         conflictsC = 0;
  vec<Lit>    learnt_clause;
  int nblevels=0;
  starts++;
  bool first = true;
  vec<Lit>    lastAssertive,intersect;
  int lastDBL=0;

  lastExtended.clear();
  for (;;){
    Clause* confl = propagate();
    if (confl != NULL){
      // CONFLICT;
      conflicts++; conflictsC++;cons++;
      if (decisionLevel() == 0) return l_False;

      first = false;

      learnt_clause.clear();
      analyze(confl, learnt_clause, backtrack_level,nblevels);


      conf4Stats++;
      nbDecisionLevelHistory.push(nblevels);
      totalSumOfDecisionLevel += nblevels;

      cancelUntil(backtrack_level);
      assert(value(learnt_clause[0]) == l_Undef);

      if (learnt_clause.size() == 1){
        uncheckedEnqueue(learnt_clause[0]);
        nbUn++; // stats
        lastAssertive.clear();
        lastExtended.clear();
      }else{
        Clause* c = Clause_new(learnt_clause, true);
        learnts.push(c);
        c->setLBD(nblevels);
        if(nblevels<=2) nbDL2++; // stats
        if(c->size()==2) nbBin++; // stats
        attachClause(*c);
        claBumpActivity(*c);
        uncheckedEnqueue(learnt_clause[0], c);

	// ALBERT
        // if(!tryToExtend(lastAssertive,lastDBL,learnt_clause,nblevels,intersect)) {
        //   lastExtended.clear();
        // }

        //        tryToExtend(lastAssertive,lastDBL,learnt_clause,nblevels,intersect);

        learnt_clause.copyTo(lastAssertive);
        lastDBL = nblevels;


      }
      varDecayActivity();
      claDecayActivity();
    }else{
      // Our dynamic restart, see the SAT09 competition compagnion paper
      if (
          ( nbDecisionLevelHistory.isvalid() &&
            ((nbDecisionLevelHistory.getavg()*0.7) > (totalSumOfDecisionLevel / conf4Stats)))) {
        nbDecisionLevelHistory.fastclear();
        progress_estimate = progressEstimate();
        cancelUntil(0);
        return l_Undef; }


      // Simplify the set of problem clauses:
      if (decisionLevel() == 0) {
        if( schedule_reduce_extend )
          reduceExtended();
        if(!simplify())
          return l_False;
      }
      //

      Lit next = lit_Undef;

      // Perform clause database reduction !
      if(cons-curRestart* nbclausesbeforereduce>=0)
        {
          curRestart = (conflicts/ nbclausesbeforereduce)+1;
          reduceDB();
          nbclausesbeforereduce += 500;
        }


      if (next == lit_Undef){
        // New variable decision:
        decisions++;
        next = pickBranchLit(polarity_mode, random_var_freq);

        if (next == lit_Undef)
          // Model found:
          return l_True;
      }

      // Increase decision level and enqueue 'next'
      assert(value(next) == l_Undef);
      newDecisionLevel();
      uncheckedEnqueue(next);
    }
  }
}


double Solver::progressEstimate() const
{
  double  progress = 0;
  double  F = 1.0 / nVars();

  for (int i = 0; i <= decisionLevel(); i++){
    int beg = i == 0 ? 0 : trail_lim[i - 1];
    int end = i == decisionLevel() ? trail.size() : trail_lim[i];
    progress += pow(F, i) * (end - beg);
  }

  return progress / nVars();
}


bool Solver::solve(const vec<Lit>& assumps)
{
  model.clear();
  conflict.clear();
  orignvars = nVars();
  schedule_reduce_extend = false;


  nbDecisionLevelHistory.initSize(100);
  totalSumOfDecisionLevel = 0;

  if (!ok) return false;

  assumps.copyTo(assumptions);

  double nof_learnts   = nClauses() * learntsize_factor;


  // Compute the first call to database clause reduction
  if(nof_learnts <nbclausesbeforereduce) {
    nbclausesbeforereduce = (nof_learnts/2 < 5000) ? 5000 : nof_learnts/2;
  }
  lbool   status        = l_Undef;

  if (verbosity >= 0){
    reportf("============================[ Search Statistics ]==============================\n");
    reportf("| Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
    reportf("|           |    Vars  Clauses extended |  nextRed  Clauses Lit/Cl |          |\n");
    reportf("===============================================================================\n");
  }

  // Search:
  while (status == l_Undef){
    if (verbosity >= 1 || starts %10==1)
      reportf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %6.3f %% |\n", (int)conflicts, order_heap.size(), nClauses(), (int)extensions, (int)nbclausesbeforereduce, nLearnts(), (double)learnts_literals/nLearnts(), progress_estimate*100), fflush(stdout);
    status = search();
  }

  if (verbosity >= 0)
    reportf("===============================================================================\n");


  if (status == l_True){
    // Extend & copy model:
    model.growTo(nVars());
    for (int i = 0; i < nVars(); i++) model[i] = value(i);
#ifndef NDEBUG
    verifyModel();
#endif
  }else{
    assert(status == l_False);
    if (conflict.size() == 0)
      ok = false;
  }

  cancelUntil(0);
  return status == l_True;
}

//=================================================================================================
// Debug methods:


void Solver::verifyModel()
{
  bool failed = false;
  for (int i = 0; i < clauses.size(); i++){
    assert(clauses[i]->mark() == 0);
    Clause& c = *clauses[i];
    for (int j = 0; j < c.size(); j++)
      if (modelValue(c[j]) == l_True)
        goto next;

    reportf("unsatisfied clause: ");
    printClause(*clauses[i]);
    reportf("\n");
    failed = true;
  next:;
  }

  assert(!failed);

  reportf("Verified %d original clauses.\n", clauses.size());
}


void Solver::checkLiteralCount()
{
  // Check that sizes are calculated correctly:
  int cnt = 0;
  for (int i = 0; i < clauses.size(); i++)
    if (clauses[i]->mark() == 0)
      cnt += clauses[i]->size();

  if ((int)clauses_literals != cnt){
    fprintf(stderr, "literal count: %d, real value = %d\n", (int)clauses_literals, cnt);
    assert((int)clauses_literals == cnt);
  }
}

