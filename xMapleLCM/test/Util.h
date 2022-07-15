#ifndef Minisat_TestUtil_h
#define Minisat_TestUtil_h

// #pragma once

#include <tr1/unordered_set>
#include <initializer_list>
#include <core/SolverTypes.h>
#include <mtl/Vec.h>

namespace Minisat {

std::tr1::unordered_set<Lit> mkLitSet(const std::initializer_list<int>& elements);
void setLitVec(vec<Lit>& v, const std::initializer_list<int>& elements);
bool requireVecEqual(Minisat::vec<Lit>& actual, Minisat::vec<Lit>& expect);
bool requireVecPrefix(Minisat::vec<Lit>& actual, Minisat::vec<Lit>& prefix);
bool requireClauseEqual(const Clause& actual, const std::initializer_list<Lit>& elements);

}

#endif