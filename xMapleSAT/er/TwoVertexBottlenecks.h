// *********************************************************
// TwoVertexBottlenecks.h -- Version 1.0, January 13, 2023
//     Author: Sam Buss (sbuss@ucsd.edu)
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
//         In this case, it should hold that predIndex.size is equal to N.
//  
// Output:
//    Return code:  0 if 3 vertex connected and thus there are no 1- or 2-vertex bottlenecks.
//                  q if there are q sets of pairs of vertices returned.
//                 -1 if not 2-vertex connected.
//                 -2 if error condition, invalid input (generates an assert)
//                
//    int NumGroups - number of groups of pairs.
//    int PairGroupInfo[]  - array of length 3*NumGroups
//              PairGroupInfo[3*k] - number of "left" elements in pairs.
//              PairGroupInfo[3*k+1] - number of "right" elements in pairs.
//              PairGroupInfo[3*k+2] - Index into GroupMembersLeft[] array.
//              PairGroupInfo[3*k+3] - Index into GroupMembersRight[] array.
//    int GroupMembersLeft[] - Holds the left vertices of vertex pairs.
//    int GroupMembersRight[] - Holds the right vertices of vertex pairs.
// All output values are stored in a TwoVertexBottlenecks class
//              
// *********************************************************

#pragma once
#ifndef TWO_VERTEX_BOTTLENECKS_H

#include <assert.h>
#include <limits>
#include <vector>

class TwoVertexBottlenecks
{
public:
    TwoVertexBottlenecks() {};

public:

    int CalcBottlenecks( int N, const int predecessors[], const int predIndex[], int maxNumGroups = INT_MAX);
    int CalcBottlenecks(const std::vector<int> predecessors, const std::vector<int> predIndex, int maxNumGroups = INT_MAX);

public:

    int NumPairGroups() const;

    int NumGroupMembersLeft(int g) const;       // Number of left members in the g-th group
    int NumGroupMembersRight(int g) const;      // Number of right members in the g-th group

    int GroupMemberLeft(int g, int k) const;    // Returns k-th left member of g-th group 
    int GroupMemberRight(int g, int k) const;   // Returns k-th right member of g-th group 
    
    // The "last" pair in a group is the furthest from the sink.
    int GroupMemberLastLeft(int g) const;       // Returns last left member of g-th group
    int GroupMemberLastRight(int g) const;      // Returns last right member of g-th group

private:
    std::vector<int> PairGroupsInfo;
    std::vector<int> GroupMembersLeft;
    std::vector<int> GroupMembersRight;

    void Clear() { PairGroupsInfo.clear(); GroupMembersLeft.clear(); GroupMembersRight.clear(); }
};

// ****************
// Returns the number of groups of pairs of 2-vertex bottlenecks.
// ***************
inline int TwoVertexBottlenecks::NumPairGroups() const {
    return PairGroupsInfo.size()/4;
}

// ******************
// Group #0 is closest to the sink.  Increasing group numbers are further from the sink.
// The "last" pair in a group is the furthest from the sink.
// ******************

// Number of left members in the g-th group
inline int TwoVertexBottlenecks::NumGroupMembersLeft(int g) const {
    assert(g < NumPairGroups());
    return PairGroupsInfo[4 * g];
}

// Number of right members in the g-th group
inline int TwoVertexBottlenecks::NumGroupMembersRight(int g) const {
    assert(g < NumPairGroups());
    return PairGroupsInfo[4 * g + 1];
}

// Returns k-th left member of g-th group
inline int TwoVertexBottlenecks::GroupMemberLeft(int g, int k) const {
    assert(g < NumPairGroups() && k < NumGroupMembersLeft(g));
    return GroupMembersLeft[PairGroupsInfo[4 * g + 2] + k];
}

// Returns k-th right member of g-th group
inline int TwoVertexBottlenecks::GroupMemberRight(int g, int k) const {
    assert(g < NumPairGroups() && k < NumGroupMembersRight(g));
    return GroupMembersRight[PairGroupsInfo[4 * g + 3] + k];
}

// Returns last left member of g-th group
inline int TwoVertexBottlenecks::GroupMemberLastLeft(int g) const {
    return GroupMemberLeft(g, NumGroupMembersLeft(g) - 1);
}

// Returns last right member of g-th group
inline int TwoVertexBottlenecks::GroupMemberLastRight(int g) const {
    return GroupMemberRight(g, NumGroupMembersRight(g) - 1);
}

// Main routine for calculating the two-vertex bottlenecks
inline int TwoVertexBottlenecks::CalcBottlenecks(
    const std::vector<int> predecessors, const std::vector<int> predIndex, int maxNumGroups)
{
    assert(predecessors.size() == predIndex.back());
    return CalcBottlenecks(predIndex.size() - 1, predecessors.data(), predIndex.data(), maxNumGroups);
}

#endif // TWO_VERTEX_BOTTLENECKS_H