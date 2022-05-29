/******************************************************************************************[Heap.h]
Copyright (c) 2022, Jonathan Chung

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

#include <cstddef>
#include <tr1/unordered_map> // Hashmap
#include <tr1/unordered_set> // Hashset
#include <utility>
#include <mtl/Vec.h>

// Template specializations for hashing
namespace std { namespace tr1 {
    // template<>
    // std::size_t std::tr1::hash<std::pair<Minisat::Lit, Minisat::Lit> >::operator()(std::pair<Minisat::Lit, Minisat::Lit> p) const {
    //     return std::size_t(p.first.x) << 32 | p.second.x;
    // }

    // template<>
    // std::size_t std::tr1::hash<Minisat::Lit>::operator()(Minisat::Lit p) const {
    //     return p.x;
    // }

    template<>
    std::size_t std::tr1::hash<std::pair<int, int> >::operator()(std::pair<int, int> p) const {
        return std::size_t(p.first) << 32 | p.second;
    }
}}

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

    static inline P mkLitPair(L a, L b) {
        return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
    }
    
public:
    inline typename PLMap::iterator find(L a, L b) { return pl_map.find(mkLitPair(a, b)); }
    inline typename PLMap::iterator end() { return pl_map.end(); }

    // Return the number of stored extension variable definitions
    inline unsigned int size() { return pl_map.size(); }

    // Return the number of extension variable definitions in which the literal participates
    inline unsigned int degree(L a) { auto it = rc_map.find(a); return it == rc_map.end() ? 0 : it->second; }

    // Check whether a pair of literals has a corresponding extension variable definition
    inline bool contains (L a, L b) { return pl_map.find(mkLitPair(a, b)) != pl_map.end(); }

    // Insert a definition for an extension variable
    // x <=> (a v b)
    void insert (L x, L a, L b) {
        const P key = mkLitPair(a, b);
        if (pl_map.find(key) != pl_map.end()) return;

        lp_map.insert(std::make_pair(x, key));
        pl_map.insert(std::make_pair(key, x));

        // Increment count for Lit a
        typename RCMap::iterator it1 = rc_map.find(a);
        if (it1 == rc_map.end()) rc_map.insert(std::make_pair(a, 1));
        else                   it1->second++;

        // Increment count for Lit b
        typename RCMap::iterator it2 = rc_map.find(b);
        if (it2 == rc_map.end()) rc_map.insert(std::make_pair(b, 1));
        else                   it2->second++;
    }

    // Delete the extension variable corresponding to a basis-literal pair
    void erase(L a, L b) {
        const P key = mkLitPair(a, b);

        // Check if the pair is in the map
        typename PLMap::iterator it = pl_map.find(key);
        if (it == pl_map.end()) return;

        // Erase from the reverse map
        lp_map.erase(it->second);

        // Decrement count for Lit a
        typename RCMap::iterator it1 = rc_map.find(a);
        if (it1->second == 1) rc_map.erase(it1);
        else                  it1->second--;

        // Decrement count for Lit b
        typename RCMap::iterator it2 = rc_map.find(b);
        if (it2->second == 1) rc_map.erase(it2);
        else                  it2->second--;

        // Erase from the forward map
        pl_map.erase(it);
    }

    // Delete a set of extension variables
    void erase(const LSet& defsToDelete) {
        for (typename LSet::const_iterator i = defsToDelete.begin(); i != defsToDelete.end(); i++) {
            typename LPMap::iterator it = lp_map.find(*i);
            if (it == lp_map.end()) continue;
            P& def = it->second;

            // Erase from the forward map
            pl_map.erase(def);

            // Decrement count for Lit a
            L a = def.first;
            typename RCMap::iterator it1 = rc_map.find(a);
            if (it1->second == 1) rc_map.erase(it1);
            else                  it1->second--;

            // Decrement count for Lit b
            L b = def.second;
            typename RCMap::iterator it2 = rc_map.find(b);
            if (it2->second == 1) rc_map.erase(it2);
            else                  it2->second--;

            // Erase from the reverse map
            lp_map.erase(it);
        }
    }

    void clear(void) {
        rc_map.clear();
        pl_map.clear();
        lp_map.clear();
    }

    void substitute(vec<L>& clause) {
        // Get indices of all basis literals (in increasing order)
        vec<int> defLitIndex;
        vec<bool> validIndex; // True if corresponding literal is in clause
        for (int i = 0; i < clause.size(); i++) {
            validIndex.push(true);
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
                if (!validIndex[defLitIndex[i]]) continue;

                // Check whether any extension variables are defined over this literal pair
                typename PLMap::iterator it = pl_map.find(mkLitPair(clause[defLitIndex[i]], clause[defLitIndex[j]]));
                if (it == end()) continue;

                // Replace the first literal with the extension literal and mark the second literal as invalid
                clause[defLitIndex[i]] = it->second;
                validIndex[defLitIndex[j]] = false;
                numReplaced++;
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