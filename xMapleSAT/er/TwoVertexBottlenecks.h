// *********************************************************
// TwoVertexBottleneck
//
// Input: A DAG with a sink and source. 
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
//    int predIndex[] - Indexed in to the predecessors array.
//          For i<N-1, if p = predIndex[i], then 
//            predecessors[p+k] is the k-th predecessor of i.
//          predIndex[N-1] is the number of entries in predecessors.
//    The DAG ordering respects the ordering of vertices, that is,
//           If j is a predecessor of i, then i<j.
//    int maxNumGroups - maximum number of pairs of vertices to return.
//           Pairs of vertices are returned in order of closest to sink first,
//           progressing towards the source.
// 
//    There is an alternate interface that uses std:vector<int> objects
//         instead of the integer arrays predecessors[] and predIndex[].  
//         In this case, it should hold that predIndex.size is equal to N.
//  
// Output:
//    Return code:  -1 if 3-vertex connected.
//                  0 if not 2-vertex connected.
//                  q if there are q sets of pairs of vertices returned.
//    int NumGroups - number of groups of pairs.
//    int PairGroupInfo[]  - array of length 3*NumGroups
//              PairGroupInfo[3*k] - number of "left" elements in pairs.
//              PairGroupInfo[3*k+1] - number of "right" elements in pairs.
//              PairGroupInfo[3*k+2] - Index into PairMembers[] array.
//    int GroupMembersArray[] - Holds the vertices in each group of vertex pairs.
// All output values are stored in a TwoVertexBottlenecks class
//              
// *********************************************************

#pragma once
#ifndef TWO_VERTEX_BOTTLENECKS_H

#include <assert.h>

class TwoVertexBottlenecks
{
public:
    TwoVertexBottlenecks() {};

    int NumPairGroups() const;

    int NumGroupMembersLeft(int g);             // Number of left members in the g-th group
    int NumGroupMembersRight(int g);            // Number of right members in the g-th group

    int GroupMemberLeft(int g, int k) const;    // Returns k-th left member of g-th group 
    int GroupMemberRight(int g, int k) const;   // Returns k-th right member of g-th group 
    
    // The "last" pair in a group is the furthest from the sink.
    int GroupMemberLastLeft(int g) const;       // Returns last left member of g-th group
    int GroupMemberLastRight(int g) const;      // Returns last right member of g-th group

private:
    std::vector<int> PairGroupsInfo;
    std::vector<int> GroupMembers;
};

// ****************
// Returns the number of groups of pairs of 2-vertex bottlenecks.
// ***************
int NumPairGroups() const {
    return PairGroupsInfo.size();
}

// ******************
// Group #0 is closest to the sink.  Increasing group numbers are further from the sink.
// The "last" pair in a group is the furthest from the sink.
// ******************

// Number of left members in the g-th group
int TwoBottlenecks::NumGroupMembersLeft(int g) const {
    assert(g < NumPairGroups());
    return PairGroupInfo[3 * g];
}

// Number of right members in the g-th group
int TwoBottlenecks::NumGroupMembersRight(int g) const {
    assert(g < NumPairGroups());
    return PairGroupInfo[3 * g + 1];
}

// Returns k-th left member of g-th group
int TwoBottlenecks::GroupMemberLeft(int g, int k) const {
    assert(g < NumPairGroups() && k < NumGroupMembersLeft(g));
    return GroupMembers[PairGroupInfo[3 * g + 2] + k];
}

// Returns k-th right member of g-th group
int TwoBottlenecks::GroupMemberRight(int g, int k) const {
    assert(g < NumPairGroups() && k < NumGroupMembersRight(g));
    return GroupMembers[PairGroupInfo[3 * g + 2] + PairGroupInfo[3 * g] + k];
}

// Returns last left member of g-th group
int TwoBottlenecks::GroupMemberLastLeft(int g) const {
    return GroupMemberLeft(NumGroupMembersLeft() - 1);
}

// Returns last rightt member of g-th group
int TwoBottlenecks::GroupMemberLastLeft(int g) const {
    return GroupMemberRight(NumGroupMembersRight() - 1);
}



#endif // TWO_VERTEX_BOTTLENECKS_H