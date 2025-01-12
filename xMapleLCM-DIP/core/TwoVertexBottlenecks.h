// *********************************************************
// TwoVertexBottlenecks.h -- 
//     Version 1.0, January 13, 2023
//     Version 2.0. October 25, 2023. Major revision, bug fix, and interface change.
// 
//     Author: Sam Buss (sbuss@ucsd.edu) with Jonathan Chung, Vijay Ganesh and Albert Oliveras
//     This code has no warranties of correctness or appropriateness.
//     May be used freely, however acknowledgement of use
//         is expected and appreciated.
//
// Input: A topologically sorted DAG with a sink and a source. 
// 
// Definition: The sink and source are k-vertex-connected
//      if there are k (directed) paths from the source to the sink which
//      are vertex-disjoint.
//
// Output: 
//   - If the sink and source are 2-vertex-connected but not 3-vertex-connected,
//     an enumeration of all distinct pairs of vertices x, y
//     such that all directed paths from the sink to the source pass through
//     at least one of x and y.
//   - Or returns that the sink is not 2-vertex-connected to the source.
//   - Or returns that the sink and the source are 3-vertex-connected.
//     
// Input contains:
//    int N - the number of vertices in the DAG.
//            Vertex 0 is the sink.  Vertex N-1 is the source.
//    int predecessors[]  - Encodes predecessors of vertices.
//    int predIndex[] - Index in to the predecessors array.
//          For i<N-1, if p = predIndex[i], then 
//            predecessors[p+k] is the k-th predecessor of i.
//          (predIndex[i+1] - predIndex[i]) is the number of predecessors of i.
//          Required: predIndex[N-1] is the number of entries in predecessors[].
//          N-1 is the source, so N-1 has no predecessors.
//    The DAG ordering respects the ordering of vertices, that is,
//           If j is a predecessor of i, then i<j.  In other words,
//           if is topologically sorted.
//    int maxNumGroups - maximum number of groups of pairs of vertices to return.
//           Pairs of vertices are returned in order of closest to sink first,
//           progressing towards the source.
// 
//    There is an alternate interface that uses std:vector<int> objects
//         instead of the integer arrays predecessors[] and predIndex[].  
//         In this case, predIndex.size must be equal to N.
//  
// Output:
//    Return code:  1 if there is at least one two vertex bottleneck.
//                  0 if 3 vertex connected and thus there are no 1- or 2-vertex bottlenecks.
//                 -1 if not 2-vertex connected (Data in SingleVertBottleneck).
//                 -2 if error condition, invalid input (generates an assert)
// 
//    VertListA - list of vertices that can be the "right" member of
//                 two vertex bottleneck.
//    VertListB - list of vertices that can be the "left" member of
//                 two vertex bottleneck.
//    Two vertex bottlenecks consist of a member of VertListA and a
//           a member of VertListB.
//    The entries in VertListA are *increasing* (going from sink to source).
//    The same holds for VertListB.
//    Each VertListA/B entry is a class VertPairInfo object:
//          - int vertNum  - vertex number in the range [1, N-1].
//          - int minPair  - a lower bound on vertices from the other VertList
//                              which can be paired with vertNum to make a bottleneck
//          - int maxPair  - an upper bound on vertices from the other VertList
//                              which can be paired with vertNum to make a bottleneck
//          - int minAncestor - any vertex in the other VertList which is greater
//                              this value is an ancestor of vertNum.
//    Note that minPair and maxPair may not actually appear in the other VertList;
//          they serve merely as lower and upper bounds. (!)
//    The same holds for minAncestor.
//    Any values vertNum1 and vertNum2 will form a two vertex bottleneck iff
//          they appear in VertListA and VertlistB and obey the indicated
//          lower and upper bounds.
//    The VertListA and VertListB values are consistent in terms of what
//          vertices form a two vertex bottleneck.  
//          That is, the VertListA minPair and maxPair values allow a pairing 
//          iff the VertListB values allow the pairing.
//    The minAncestor information is not redundant however.
// 
//    GetPathA() and GetPathB():
//       PathA & PathB - returned as complete vertex-disjoint complete paths 
//          from the vertex 0 to vertex N-1.
//          All TVB pairs have one member from PathA and one member
//          from PathB.        
// *********************************************************

#pragma once
#ifndef TWO_VERTEX_BOTTLENECKS_H

#include <assert.h>
#include <limits>
#include <climits>
#include <vector>

class TwoVertexBottlenecks
{
public:
    TwoVertexBottlenecks() {};

public:

    int CalcBottlenecks( int N, const int predecessors[], const int predIndex[], int maxNumGroups = INT_MAX);
    int CalcBottlenecks(const std::vector<int> predecessors, const std::vector<int> predIndex, int maxNumGroups = INT_MAX);

public:
    class VertPairInfo {
    public:
        VertPairInfo(int vert, int lowerBdPair, int upperBdPair, int minAncestorOther)
            : vertNum(vert), minPairIdx(lowerBdPair), maxPairIdx(upperBdPair), minAncestor(minAncestorOther) {}

        int vertNum;        // The index of a vertex as potential member of 2-vertex bottleneck
        int minPairIdx;     // Index of the first vertex in the other path whigh forms a TVB with this vertex.
        int maxPairIdx;     // Index of the last vertex in the other path whigh forms a TVB with this vertex.
        int minAncestor;    // Lower bound on verts in the other VertPairInfo array which are ancestors
    };

public:
    const std::vector<VertPairInfo>& GetVertListA() const { return VertListA; }
    const std::vector<VertPairInfo>& GetVertListB() const { return VertListB; }

    const std::vector<int>& GetPathA() const { return PathA; }
    const std::vector<int>& GetPathB() const { return PathB; }

    int LengthListA() const { return (int)VertListA.size(); }
    int LengthListB() const { return (int)VertListA.size(); }
    const VertPairInfo& GetListA(int n) const { return VertListA[n]; }
    const VertPairInfo& GetListB(int n) const { return VertListB[n]; }

public:    // If not 2-connected (return code -1), this is the single vertex bottleneck closest to the sink
    int SingleVertBottleneck;

private:
    // VertListA and VertListB hold the two lists of members of PathA and PathB
    //   that are part of a DIPs. They also old the minVert and maxVert info
    //   characterizing which vertices they form pairs with.
    std::vector<VertPairInfo> VertListA;
    std::vector<VertPairInfo> VertListB;

    // PathA and PathB are the complete vertex-disjoint paths.
    //   They also hold the start and end vertices N-1 and 0.
    // PathA and PathB are reporting information only. Many applications
    //   may not need them.  In any event, the complete paths are not
    //   unique. Nor are they needed internally for computation by CalcBottlenecks.
    std::vector<int> PathA;
    std::vector<int> PathB;

    void Clear() { VertListA.clear(); VertListB.clear(); }

};


// Main routine for calculating the two-vertex bottlenecks
//   This is just a wrapper to accept std::vector's as inputs.
inline int TwoVertexBottlenecks::CalcBottlenecks(
    const std::vector<int> predecessors, const std::vector<int> predIndex, int maxNumGroups)
{
    assert((int)predecessors.size() == predIndex.back());
    return CalcBottlenecks((int)predIndex.size(), predecessors.data(), predIndex.data(), maxNumGroups);
}

#endif // TWO_VERTEX_BOTTLENECKS_H