#ifndef Minisat_TestUtil_h
#define Minisat_TestUtil_h

// #pragma once

#include <initializer_list>
#include <tr1/unordered_set>
#include <core/SolverTypes.h>
#include <mtl/Vec.h>

namespace Minisat {

// Utility methods
std::tr1::unordered_set<Lit> mkLitSet(const std::initializer_list<int>& elements);
void setLitVec(vec<Lit>& v, const std::initializer_list<int>& elements);
void clause2Vec(vec<Lit>& v, const Clause& c);
template <typename T> void setVec(vec<T>& v, const std::initializer_list<T>& elements);
template <typename T> void setVec(vec<T>& v, const std::tr1::unordered_set<T>& c);

// Matchers
class VecPrefix;
class VecEqual;
class VecEqualUnordered;

// Vector prefix matcher and builder function
VecPrefix vecPrefix(const Minisat::vec<Lit>& prefix);
class VecPrefix : public Catch::MatcherBase<Minisat::vec<Lit>> {
public:
    VecPrefix(const Minisat::vec<Lit>& prefix) : m_prefix(prefix) {}
    virtual bool match(const Minisat::vec<Lit>& actual) const override;
    virtual std::string describe() const override;
private:
    const Minisat::vec<Lit>& m_prefix;
};

// Vector equality matcher and builder function
VecEqual vecEqual(const Minisat::vec<Lit>& expect);
class VecEqual : public Catch::MatcherBase<Minisat::vec<Lit>> {
public:
    VecEqual(const Minisat::vec<Lit>& expect) : m_expect(expect) {}
    virtual bool match(const Minisat::vec<Lit>& actual) const override;
    virtual std::string describe() const override;
private:
    const Minisat::vec<Lit>& m_expect;
};

// Vector unordered equality matcher and builder function
VecEqualUnordered vecEqualUnordered(const Minisat::vec<Lit>& expect);
class VecEqualUnordered : public Catch::MatcherBase<Minisat::vec<Lit>> {
public:
    VecEqualUnordered(const Minisat::vec<Lit>& expect) : m_expect(expect) {}
    virtual bool match(const Minisat::vec<Lit>& actual) const override;
    virtual std::string describe() const override;
private:
    const Minisat::vec<Lit>& m_expect;
};

}

#endif