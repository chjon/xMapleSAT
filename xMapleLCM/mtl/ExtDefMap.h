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
    RCMap lits;
    
    // Map of definitions to extension variable
    PLMap map1;
    // Map of extension variable to definition
    LPMap map2;

    static inline P mkLitPair(L a, L b) {
        return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
    }
    
public:
    inline typename PLMap::iterator find(L a, L b) { return map1.find(mkLitPair(a, b)); }
    inline typename PLMap::iterator end() { return map1.end(); }

    // Return the number of stored extension variable definitions
    inline unsigned int size() { return map1.size(); }

    // Return the number of extension variable definitions in which the literal participates
    inline unsigned int degree(L a) { auto it = lits.find(a); return it == lits.end() ? 0 : it->second; }

    inline bool contains (L a     ) { return lits.find(a)               != lits.end(); }
    inline bool contains (L a, L b) { return map1.find(mkLitPair(a, b)) != map1.end(); }

    void insert (L x, L a, L b) {
        const P key = mkLitPair(a, b);
        if (map1.find(key) != map1.end()) return;

        map2.insert(std::make_pair(x, key));
        map1.insert(std::make_pair(key, x));

        // Increment count for Lit a
        typename RCMap::iterator it1 = lits.find(a);
        if (it1 == lits.end()) lits.insert(std::make_pair(a, 1));
        else                   it1->second++;

        // Increment count for Lit b
        typename RCMap::iterator it2 = lits.find(b);
        if (it2 == lits.end()) lits.insert(std::make_pair(b, 1));
        else                   it2->second++;
    }

    void erase(L a, L b) {
        const P key = mkLitPair(a, b);

        // Check if the pair is in the map
        typename PLMap::iterator it = map1.find(key);
        if (it == map1.end()) return;

        // Erase from the reverse map
        map2.erase(it->second);

        // Decrement count for Lit a
        typename RCMap::iterator it1 = lits.find(a);
        if (it1->second == 1) lits.erase(it1);
        else                  it1->second--;

        // Decrement count for Lit b
        typename RCMap::iterator it2 = lits.find(b);
        if (it2->second == 1) lits.erase(it2);
        else                  it2->second--;

        // Erase from the forward map
        map1.erase(it);
    }

    void erase(const LSet& defsToDelete) {
        for (typename LSet::const_iterator i = defsToDelete.begin(); i != defsToDelete.end(); i++) {
            typename LPMap::iterator it = map2.find(*i);
            if (it == map2.end()) continue;
            P& def = it->second;

            // Erase from the forward map
            map1.erase(def);

            // Decrement count for Lit a
            L a = def.first;
            typename RCMap::iterator it1 = lits.find(a);
            if (it1->second == 1) lits.erase(it1);
            else                  it1->second--;

            // Decrement count for Lit b
            L b = def.second;
            typename RCMap::iterator it2 = lits.find(b);
            if (it2->second == 1) lits.erase(it2);
            else                  it2->second--;

            // Erase from the reverse map
            map2.erase(it);
        }
    }

    void clear(void) {
        lits.clear();
        map1.clear();
        map2.clear();
    }
};

}

#endif