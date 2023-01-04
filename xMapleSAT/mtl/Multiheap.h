#ifndef Minisat_Multiheap_h
#define Minisat_Multiheap_h

#include "mtl/Heap.h"
#include "mtl/Vec.h"

namespace Minisat {
template<class LT>
class Multiheap {
    Heap<LT> m_root_heap;
    LT m_sub_lt;
    vec< Heap<LT> > m_heaps;
    const vec<uint64_t>& m_root_map;

    int m_size;
public:
    Multiheap(const LT& root_lt, const LT& sub_lt, const vec<uint64_t>& root_map)
        : m_root_heap(root_lt)
        , m_sub_lt(sub_lt)
        , m_root_map(root_map)
        , m_size(0)
    {}

    int  size      ()          const { return m_size; }
    bool empty     ()          const { return m_size == 0; }
    bool inHeap    (int n)     const {
        const int rootKey = m_root_map[n];
        if (!m_root_heap.inHeap(rootKey)) return false;
        return m_heaps[rootKey].inHeap(n);
    }
    int operator[](int index) const {
        assert(index < m_size);
        int i = 0;
        while (index >= m_heaps[i].size())
            index -= m_heaps[i++].size();
        return m_heaps[i][index];
    }

    void decrease  (int n) { const int rootKey = m_root_map[n]; m_heaps[rootKey].decrease(n); m_root_heap.decrease(rootKey); }
    void increase  (int n) { const int rootKey = m_root_map[n]; m_heaps[rootKey].increase(n); m_root_heap.increase(rootKey); }

    // Safe variant of insert/decrease/increase:
    void update(int n) {
        if (!inHeap(n))
            insert(n);
        else {
            decrease(n);
            increase(n);
        }
    }

    void insert(int n) {
        const int rootKey = m_root_map[n];
        if (!m_root_heap.inHeap(rootKey)) m_root_heap.insert(rootKey);
        while (rootKey >= m_heaps.size())
            m_heaps.push_new(reinterpret_cast<void*>(&m_sub_lt));
        m_heaps[rootKey].insert(n);
        m_size++;
    }

    int  removeMin() {
        int i = m_root_heap[0];
        int x = m_heaps[i].removeMin();
        if (m_heaps[i].size() == 0)
            m_root_heap.removeMin();
        m_size--;
        return x; 
    }

    // Rebuild the heap from scratch, using the elements in 'ns':
    void build(vec<int>& ns) {
        // Clear root heap
        m_root_heap.clear();
        for (int i = 0; i < m_heaps.size(); i++)
            m_heaps[i].clear();

        // Rebuild heap
        for (int i = 0; i < ns.size(); i++)
            insert(ns[i]);

        m_size = ns.size();
    }

    void clear(bool dealloc = false) 
    { 
        m_root_heap.clear(dealloc);
        for (int i = 0; i < m_heaps.size(); i++)
            m_heaps[i].clear(dealloc);

        m_size = 0;
    }
};
}

#endif