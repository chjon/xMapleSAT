/******************************************************************************************[ExtDefMap.h]
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

#include <stdio.h>
#include <sstream>
#include <tr1/unordered_map>

#include "er/ERManager.h"
#include "test/catch.hpp"
#include "test/Util.h"

namespace Minisat {

std::tr1::unordered_set<Lit> mkLitSet(const std::initializer_list<int>& elements) {
    std::tr1::unordered_set<Lit> s;
    for (const auto element : elements) s.insert(mkLit(element < 0 ? -element : element, element < 0));
    return s;
}

void setLitVec(vec<Lit>& v, const std::initializer_list<int>& elements) {
    v.clear();
    for (const auto element : elements) v.push(mkLit(element < 0 ? -element : element, element < 0));
}

template <typename T> void setVec(vec<T>& v, const std::initializer_list<T>& elements) {
    v.clear();
    for (const auto element : elements) v.push(element);
}
template void setVec(vec<Lit>&, const std::initializer_list<Lit>&);
template void setVec(vec<int>&, const std::initializer_list<int>&);

void clause2Vec(vec<Lit>& v, const Clause& c) {
    v.clear();
    for (int i = 0; i < c.size(); i++) v.push(c[i]);
}

template <typename T> void setVec(vec<T>& v, const std::tr1::unordered_set<T>& elements) {
    v.clear();
    for (const auto element : elements) v.push(element);
}
template void setVec(vec<Lit>&, const std::tr1::unordered_set<Lit>&);
template void setVec(vec<int>&, const std::tr1::unordered_set<int>&);

void setVariables(AssignmentTrail& at, int i_undef, int i_max, int numVars) {
    at.cancelUntil(0);
    for (int i = 0; i < numVars; i++)
        if (i != i_undef && i != i_max) {
            at.newDecisionLevel();
            at.assign(mkLit(i, true));
        }
    at.newDecisionLevel();
    at.assign(mkLit(i_max, true));
}

///////////////////////////
// Vector prefix matcher //
///////////////////////////

VecPrefix vecPrefix(const Minisat::vec<Lit>& prefix) { return VecPrefix(prefix); }

bool VecPrefix::match(const Minisat::vec<Lit>& actual) const {
    if (actual.size() < m_prefix.size()) return false;
    for (int i = 0; i < m_prefix.size(); i++)
        if (actual[i] != m_prefix[i])
            return false;
    return true;
}

std::string VecPrefix::describe() const {
    std::ostringstream ss;
    ss << "has prefix " << m_prefix;
    return ss.str();
}

/////////////////////////////
// Vector equality matcher //
/////////////////////////////

VecEqual vecEqual(const Minisat::vec<Lit>& expect) { return VecEqual(expect); }

bool VecEqual::match(const Minisat::vec<Lit>& actual) const {
    if (actual.size() != m_expect.size()) return false;
    for (int i = 0; i < actual.size(); i++)
        if (actual[i] != m_expect[i])
            return false;
    return true;
}

std::string VecEqual::describe() const {
    std::ostringstream ss;
    ss << "is equal to " << m_expect;
    return ss.str();
}

///////////////////////////////////////
// Vector unordered equality matcher //
///////////////////////////////////////

template <typename T>
VecEqualUnordered<T> vecEqualUnordered(const Minisat::vec<T>& expect) { return VecEqualUnordered<T>(expect); }
template VecEqualUnordered<Lit> vecEqualUnordered(const Minisat::vec<Lit>&);
template VecEqualUnordered<int> vecEqualUnordered(const Minisat::vec<int>&);

template <typename T>
bool VecEqualUnordered<T>::match(const Minisat::vec<T>& actual) const {
    std::tr1::unordered_map<T, int> count;
    if (actual.size() != m_expect.size()) return false;
    for (int i = 0; i < m_expect.size(); i++) {
        auto it = count.find(m_expect[i]);
        if (it == count.end()) count.insert(std::make_pair(m_expect[i], 1));
        else it->second++;
    }

    for (int i = 0; i < actual.size(); i++) {
        auto it = count.find(actual[i]);
        if (it == count.end()) return false;
        else it->second--;
    }

    for (auto i = count.begin(); i != count.end(); i++)
        if (i->second != 0) return false;
    return true;
}
template bool VecEqualUnordered<Lit>::match(const Minisat::vec<Lit>&) const;
template bool VecEqualUnordered<int>::match(const Minisat::vec<int>&) const;

template <typename T>
std::string VecEqualUnordered<T>::describe() const {
    std::ostringstream ss;
    ss << "is equal to a permutation of " << m_expect;
    return ss.str();
}
template std::string VecEqualUnordered<Lit>::describe() const;
template std::string VecEqualUnordered<int>::describe() const;

///////////////////////////////
// ExtDef uniqueness matcher //
///////////////////////////////

ExtDefUnique extDefUnique() { return ExtDefUnique(); }

bool ExtDefUnique::match(const std::vector<ExtDef>& v) const {
    std::tr1::unordered_set<Lit> seenVar;
    std::tr1::unordered_set< std::pair<Lit, Lit> > seenDef;
    for (ExtDef def : v) {
        std::pair<Lit, Lit> p = (def.a < def.b) ? std::make_pair(def.a, def.b) : std::make_pair(def.b, def.a);
        if (seenVar.find(def.x) != seenVar.end()) return false;
        if (seenDef.find(p) != seenDef.end()) return false;
        seenVar.insert(def.x); seenDef.insert(p);
    }
    return true;
}

std::string ExtDefUnique::describe() const {
    std::ostringstream ss;
    ss << "only contains unique definitions";
    return ss.str();
}

/////////////////////////////////
// Watcher correctness matcher //
/////////////////////////////////

WatchersCorrect watchersCorrect(Minisat::OccLists<Lit, vec<Watcher>, WatcherDeleted>& ws, CRef cr) { return WatchersCorrect(ws, cr); }

bool WatchersCorrect::match(const Minisat::vec<Lit>& c) const {
    m_failure = 0;

    if (m_ws[ c[0]].size() != 0) m_failure |= FailureCase::EXPECT_ZERO_WATCH0  ;
    if (m_ws[ c[1]].size() != 0) m_failure |= FailureCase::EXPECT_ZERO_WATCH1  ;
    if (m_ws[~c[0]].size() != 1) m_failure |= FailureCase::EXPECT_SINGLE_WATCH0;
    if (m_ws[~c[1]].size() != 1) m_failure |= FailureCase::EXPECT_SINGLE_WATCH0;

    for (int i = 2; i < c.size(); i++) {
        if (m_ws[ c[i]].size() != 0) m_failure |= FailureCase::EXPECT_ZERO_OTHER;
        if (m_ws[~c[i]].size() != 0) m_failure |= FailureCase::EXPECT_ZERO_OTHER;
    }

    if (!find(m_ws[~c[0]], Watcher(m_cr, c[1]))) m_failure |= FailureCase::EXPECT_FOUND_WATCH0;
    if (!find(m_ws[~c[1]], Watcher(m_cr, c[0]))) m_failure |= FailureCase::EXPECT_FOUND_WATCH1;

    return m_failure == 0;
}

std::string WatchersCorrect::describe() const {
    std::ostringstream ss;

    if (m_failure & EXPECT_ZERO_WATCH0  ) ss << "has 0 watchers for c[0]; ";
    if (m_failure & EXPECT_ZERO_WATCH1  ) ss << "has 0 watchers for c[1]; ";
    if (m_failure & EXPECT_SINGLE_WATCH0) ss << "has 1 watcher for ~c[0]; ";
    if (m_failure & EXPECT_SINGLE_WATCH1) ss << "has 1 watcher for ~c[1]; ";
    if (m_failure & EXPECT_ZERO_OTHER   ) ss << "has no watchers for literals other than c[0] and c[1]; ";
    if (m_failure & EXPECT_FOUND_WATCH0 ) ss << "~c[0] watches c[1]; ";
    if (m_failure & EXPECT_FOUND_WATCH0 ) ss << "~c[1] watches c[0]; ";
    
    return ss.str();
}

}