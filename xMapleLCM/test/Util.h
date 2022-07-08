#ifndef Minisat_TestUtil_h
#define Minisat_TestUtil_h

// #pragma once

#include <tr1/unordered_set>
#include <initializer_list>
#include <core/SolverTypes.h>
#include <mtl/Vec.h>

namespace Minisat {

std::tr1::unordered_set<Lit> mkLitSet(std::initializer_list<int> elements);
void setLitVec(vec<Lit>& v, std::initializer_list<int> elements);
bool requireVecEqual(Minisat::vec<Lit>& actual, Minisat::vec<Lit>& expect);
bool requireVecPrefix(Minisat::vec<Lit>& actual, Minisat::vec<Lit>& prefix);

}

#endif