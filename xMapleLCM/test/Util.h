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

#ifndef Minisat_TestUtil_h
#define Minisat_TestUtil_h

// #pragma once

#include <initializer_list>
#include <tr1/unordered_set>
#include <core/Solver.h>
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
class ExtDefUnique;
class WatchersCorrect;

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

// ExtDef uniqueness matcher and builder function
ExtDefUnique extDefUnique();
class ExtDefUnique : public Catch::MatcherBase<const std::vector<ExtDef>> {
public:
    virtual bool match(const std::vector<ExtDef>&) const override;
    virtual std::string describe() const override;
};

// Watcher correctness matcher
WatchersCorrect watchersCorrect(Minisat::OccLists<Lit, vec<Solver::Watcher>, Solver::WatcherDeleted>& ws, CRef cr);
class WatchersCorrect : public Catch::MatcherBase<Minisat::vec<Lit>> {
public:
    WatchersCorrect(Minisat::OccLists<Lit, vec<Solver::Watcher>, Solver::WatcherDeleted>& ws, CRef cr) : m_ws(ws), m_cr(cr) {}
    virtual bool match(const Minisat::vec<Lit>&) const override;
    virtual std::string describe() const override;
private:
    enum FailureCase : int {
        EXPECT_ZERO_WATCH0   = 1 << 0,
        EXPECT_ZERO_WATCH1   = 1 << 1,
        EXPECT_SINGLE_WATCH0 = 1 << 2,
        EXPECT_SINGLE_WATCH1 = 1 << 3,
        EXPECT_FOUND_WATCH0  = 1 << 4,
        EXPECT_FOUND_WATCH1  = 1 << 6,
        EXPECT_ZERO_OTHER    = 1 << 7,
    };

    Minisat::OccLists<Lit, vec<Solver::Watcher>, Solver::WatcherDeleted>& m_ws;
    const CRef m_cr;
    mutable int m_failure = 0;
};

}

#endif