#include <stdio.h>
#include <sstream>
#include <tr1/unordered_map>
#include "catch.hpp"
#include "Util.h"

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

void clause2Vec(vec<Lit>& v, const Clause& c) {
    v.clear();
    for (int i = 0; i < c.size(); i++) v.push(c[i]);
}

template <typename T> void setVec(vec<T>& v, const std::tr1::unordered_set<T>& elements) {
    v.clear();
    for (const auto element : elements) v.push(element);
}
template void setVec(vec<Lit>&, const std::tr1::unordered_set<Lit>&);

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

VecEqualUnordered vecEqualUnordered(const Minisat::vec<Lit>& expect) { return VecEqualUnordered(expect); }

bool VecEqualUnordered::match(const Minisat::vec<Lit>& actual) const {
    std::tr1::unordered_map<Lit, int> count;
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

std::string VecEqualUnordered::describe() const {
    std::ostringstream ss;
    ss << "is equal to a permutation of " << m_expect;
    return ss.str();
}

}