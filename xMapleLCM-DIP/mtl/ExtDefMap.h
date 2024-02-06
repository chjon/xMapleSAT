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

#ifndef xMaplesat_ExtDefMap_h
#define xMaplesat_ExtDefMap_h

#ifndef USE_NONBASIS_VAR_SET
#define USE_NONBASIS_VAR_SET false
#endif

#include <cstddef>
#include <tr1/unordered_map> // Hashmap
#include <tr1/unordered_set> // Hashset
#include <utility>
#include <mtl/Vec.h>

namespace Minisat {

template<class L> // Literal - must implement operator<
class ExtDefMap {
private:
    // Pair
    using P = std::pair<L, L>;

    // Pair-literal map
    using PLMap = std::tr1::unordered_map<P, L>;

    // Literal-pair map
    using LPMap = std::tr1::unordered_map<L, P>;

    // Reference-count map
    using RCMap = std::tr1::unordered_map<L, int>;

    // Literal set
    using LSet = std::tr1::unordered_set<L>;

    // Count of all literals that appear in an extension variable definition
    RCMap rc_map;
    
    // Map of definitions to extension variable
    PLMap pl_map;
    
    // Map of extension variable to definition
    LPMap lp_map;

#if USE_NONBASIS_VAR_SET
    // Extension variables which do not participate in definitions
    LitSet nonbasis_lits;
#endif

    static inline P mkLitPair(L a, L b) {
        return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
    }

    // Vectors used for substitution -- allocated here to avoid repeated memory allocation
    mutable vec<int> defLitIndex; // The indices of all the basis literals in a clause
    mutable vec<bool> validIndex; // Flag for whether a literal belongs to the final clause

public:
    inline typename LPMap::const_iterator find(L x) const { return lp_map.find(x); }
    inline typename PLMap::const_iterator find(L a, L b) const { return pl_map.find(mkLitPair(a, b)); }
    //inline typename PLMap::iterator end() { return pl_map.end(); }

    // Return the number of stored extension variable definitions
    inline unsigned int size() const { return pl_map.size(); }

    // Return the number of extension variable definitions in which the literal participates
    inline unsigned int degree(L a) const { auto it = rc_map.find(a); return it == rc_map.end() ? 0 : it->second; }

    // Check whether a pair of literals has a corresponding extension variable definition
    inline bool containsPair (L a, L b) const { return pl_map.find(mkLitPair(a, b)) != pl_map.end(); }

    // Returns whether a pair of literals has a corresponding exension variable definition and if so, which one
    std::pair<bool,Lit> findPair(L a, L b) const {
      auto it = find(a,b);
      if (it == pl_map.end()) return {false,mkLit(1,false)};
      else return {true,it->second};
    }
  
    // Check whether a variable is a current extension variable
    inline bool containsExt (L x) const { return lp_map.find(x) != lp_map.end(); }

    // Get set of non-basis extension variables
    inline LSet getNonbasisExtLits () const {
    #if USE_NONBASIS_VAR_SET
        return nonbasis_lits;
    #else
        LSet tmp;
        for (auto it = lp_map.begin(); it != lp_map.end(); it++) {
            const Lit x = it->first;
            if (degree(x) == 0 && degree(~x) == 0) tmp.insert(x);
        }
        return tmp;
    #endif
    }

    // Insert a definition for an extension variable
    // x <=> (a v b)
    void insert (L x, L a, L b) {
        const P key = mkLitPair(a, b);
        if (pl_map.find(key) != pl_map.end()) return;

        auto lp_pair = lp_map.insert(std::make_pair(x, key));
        auto pl_pair = pl_map.insert(std::make_pair(key, x));
        assert(lp_pair.second && pl_pair.second);
    #if USE_NONBASIS_VAR_SET
        nonbasis_lits.insert(x);
    #endif

        // Increment count for Lit a
        typename RCMap::iterator it1 = rc_map.find(a);
        if (it1 == rc_map.end()) {
            rc_map.insert(std::make_pair(a, 1));
        #if USE_NONBASIS_VAR_SET
            nonbasis_lits.erase(a); nonbasis_lits.erase(~a);
        #endif
        } else it1->second++;

        // Increment count for Lit b
        typename RCMap::iterator it2 = rc_map.find(b);
        if (it2 == rc_map.end()) {
            rc_map.insert(std::make_pair(b, 1));
        #if USE_NONBASIS_VAR_SET
            nonbasis_lits.erase(b); nonbasis_lits.erase(~b);
        #endif
        } else it2->second++;
    }

    // Delete the extension variable corresponding to a basis-literal pair
    void erase(L a, L b) {
        const P key = mkLitPair(a, b);

        // Check if the pair is in the map
        typename PLMap::iterator it = pl_map.find(key);
        if (it == pl_map.end()) return;
        Lit x = it->second;

    #if USE_NONBASIS_VAR_SET
        // Erase from the non-basis lit set
        assert(nonbasis_lits.find(x) != nonbasis_lits.end());
        nonbasis_lits.erase(x);
    #endif

        // Erase from the reverse map
        lp_map.erase(x);

        // Decrement count for Lit a
        typename RCMap::iterator it1 = rc_map.find(a);
        if (it1->second == 1) {
            rc_map.erase(it1);
        #if USE_NONBASIS_VAR_SET
            if (lp_map.find( a) != lp_map.end()) nonbasis_lits.insert( a);
            if (lp_map.find(~a) != lp_map.end()) nonbasis_lits.insert(~a);
        #endif
        } else it1->second--;

        // Decrement count for Lit b
        typename RCMap::iterator it2 = rc_map.find(b);
        if (it2->second == 1) {
            rc_map.erase(it2);
        #if USE_NONBASIS_VAR_SET
            if (lp_map.find( b) != lp_map.end()) nonbasis_lits.insert( b);
            if (lp_map.find(~b) != lp_map.end()) nonbasis_lits.insert(~b);
        #endif
        } else it2->second--;

        // Erase from the forward map
        pl_map.erase(it);
    }

    // Delete an extension variable
    void erase(L l) {
        typename LPMap::iterator it = lp_map.find(l);
        if (it == lp_map.end()) return;
        P& def = it->second;

    #if USE_NONBASIS_VAR_SET
        // Erase from the non-basis lit set
        assert(nonbasis_lits.find(l) != nonbasis_lits.end());
        nonbasis_lits.erase(l);
    #endif

        // Erase from the forward map
        pl_map.erase(def);

        // Decrement count for Lit a
        L a = def.first;
        typename RCMap::iterator it1 = rc_map.find(a);
        if (it1->second == 1) {
            rc_map.erase(it1);
        #if USE_NONBASIS_VAR_SET
            if (lp_map.find( a) != lp_map.end()) nonbasis_lits.insert( a);
            if (lp_map.find(~a) != lp_map.end()) nonbasis_lits.insert(~a);
        #endif
        } else it1->second--;

        // Decrement count for Lit b
        L b = def.second;
        typename RCMap::iterator it2 = rc_map.find(b);
        if (it2->second == 1) {
            rc_map.erase(it2);
        #if USE_NONBASIS_VAR_SET
            if (lp_map.find( b) != lp_map.end()) nonbasis_lits.insert( b);
            if (lp_map.find(~b) != lp_map.end()) nonbasis_lits.insert(~b);
        #endif
        } else it2->second--;

        // Erase from the reverse map
        lp_map.erase(it);
    }

    void clear(void) {
        rc_map.clear();
        pl_map.clear();
        lp_map.clear();
    #if USE_NONBASIS_VAR_SET
        nonbasis_lits.clear();
    #endif
    }

    void absorb(vec<L>& clause) const {
        std::tr1::unordered_map<L, int> basis; // Map basis literals to the index of their extension lits
        validIndex.clear();

        // Find all the literals in the definitions of extension literals in the clause
        for (int i = 0; i < clause.size(); i++) {
            validIndex.push(true);

            // Check if variable is an extension variable
            auto it = lp_map.find(clause[i]);
            if (it == lp_map.end()) {
                it = lp_map.find(~clause[i]);
                if (it == lp_map.end()) continue;
            }

            const std::pair<L, L> ab = it->second;
            basis.insert(std::make_pair(ab.first , i));
            basis.insert(std::make_pair(ab.second, i));
        }

        if (basis.size() > 0) {
            // Mark redundant literals
            for (int i = 0; i < clause.size(); i++) {
                // Delete both a and b if the clause contains x
                if (basis.find(clause[i]) != basis.end()) validIndex[i] = false;

                // Delete -x if the clause contains -a or -b
                const auto it = basis.find(~clause[i]);
                if (it != basis.end() && sign(clause[it->second])) validIndex[it->second] = false;
            }

            // Sweep through clause, removing redundant literals
            // Never remove the literal at index 0
            int i, j;
            for (i = 1, j = 1; i < clause.size(); i++) {
                clause[j] = clause[i];
                j += validIndex[i];
            }
            clause.shrink(clause.size() - j);
        }
    }

    void substitute(vec<L>& clause, vec<L>& extLits) const {
        // Get indices of all basis literals (in increasing order)
        LSet lits;
        defLitIndex.clear();
        validIndex.clear();
        for (int i = 0; i < clause.size(); i++) {
            validIndex.push(true);
            lits.insert(clause[i]);
            // Check if any extension variables are defined over this literal
            if (rc_map.find(clause[i]) != rc_map.end()) {
                defLitIndex.push(i);
            }
        }
        if (defLitIndex.size() <= 1) return;

        // Check each pair of basis literals for a corresponding extension variable
        int numReplaced = 0;
        for (int i = 0; i < defLitIndex.size(); i++) {
            if (!validIndex[defLitIndex[i]]) continue;
            for (int j = i + 1; j < defLitIndex.size(); j++) {
                if (!validIndex[defLitIndex[j]]) continue;

                // Check whether any extension variables are defined over this literal pair
                typename PLMap::const_iterator it = pl_map.find(mkLitPair(clause[defLitIndex[i]], clause[defLitIndex[j]]));
                if (it == pl_map.end()) continue;

                // Check whether the extension literal is already present in the clause
                if (lits.find(it->second) != lits.end()) continue;

                // Replace the first literal with the extension literal and mark the second literal as invalid
                clause[defLitIndex[i]] = it->second;
                validIndex[defLitIndex[j]] = false;
                numReplaced++;
                extLits.push(it->second);
                break;
            }
        }

        // Generate reduced learnt clause
        if (numReplaced > 0) {
            for (int i = 0, j = 0; i < clause.size(); i++) {
                // Impl1:
                // Naive implementation:
                // if (validIndex[i]) clause[j++] = clause[i];

                // Impl2:
                // Avoid branching
                // Implementation conducive to SIMD instructions (similar to Stream VByte)
                clause[j] = clause[i];
                j += validIndex[i];
            }
            clause.shrink(numReplaced);
        }

        // Check whether the clause contains multiple of the same variable (should never happen)
        // LSet clauseLits;
        // for (int i = 0; i < clause.size(); i++) {
        //     assert(clauseLits.find(clause[i]) == clauseLits.end());
        //     clauseLits.insert(clause[i]);
        // }
    }
};

}

#endif
