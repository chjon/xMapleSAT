
// TwoVertexBottlenecks.cpp -- Version 1.0, January 13, 2023
//     Author: Sam Buss (sbuss@ucsd.edu)
//     This code has no warranties of correctness or appropriateness.
//     May be used freely, however acknowledgement of use 
//         is expected and appreciated.
// 
//  ************************************************
//  * Main routines for TwoVectexBottlenecks       *
//  *                                              *
//  * See TwoVertexBottleneck.h for usage.         *
//  ************************************************

#include "TwoVertexBottlenecks.h"
#include <algorithm>
#include <limits>

// ********************************************
// Internal data structures.
// ********************************************

//   vertInfo-- uSed in Phases A, B, C
//   Tracks for each vertex: whether it is on a path
//           and whether it has been visited the depth first searches.

class vertInfo {
public:
    int status = 0;     // Status
    int succ = -1;      // Successor node

    bool TestOnFirst() const { return status & 0x01; }
    void SetOnFirst() { status = status | 0x01; }
    void ResetOnFirst() { status = (status & ~0x01); }

    bool TestOnSecond() const { return status & 0x02; }
    void SetOnSecond() { status = status | 0x02; }

    int TestOnPath() const { return (status & 0x03); }   // 0, 1, 2, 3 - If on neither path, path one, path two or both.
    bool TestOnPath(int leftFlag) { return status & (0x01 << leftFlag); }

    bool TestReachedPhaseB() const { return status & 0x04; }
    bool ReachedInPhaseB() {
        bool ret = status & 0x04;
        status = status | 0x04;
        return ret;
    }

    // leftFlag will be 0 or 1 for working on reachability 
    //     from left path vertices or from right path vertices (respectively)
    bool ReachedInPhaseC( int leftFlag ) {
        bool ret = status & (0x08 << leftFlag);
        status = status | (0x08 << leftFlag);
        return ret;
    }
};

class pivotInfo {
    // In the array of pathJumpingInfo objects:
    //     First PivotVert 0. (Sink)   Last PivotVert is N-1. (Source)
    //     MaxReachOtherPath is *non-decreasing*

public:

    pivotInfo(int vert, int maxReach) { PivotVert = vert; MaxReachOtherPath = maxReach; }

    int PivotVert;  
    int MaxReachOtherPath;
};

int TwoVertexBottlenecks::CalcBottlenecks(int N, const int predecessors[], const int predIndex[], int maxNumGroups)
{
    // Allocate vertex info array.
    std::vector<vertInfo> verts(N);
    
    Clear();    // Initialize output info 

    // ***************************************************
    // Phase A: Construct a path from the sink to the source (in direction of reverse edges)
    // ***************************************************

    int curVert = 0;        // Start at the sink
    verts[0].SetOnFirst();
    while (curVert != N-1) {
        int nextVert = predecessors[predIndex[curVert]];
        if (nextVert <= curVert || nextVert > N-1) {
            assert("Invalid predecessor" && false);
            return -2;
        }
        verts[nextVert].succ = curVert;
        verts[nextVert].SetOnFirst();
        curVert = nextVert;
    }

    // ****************************************************
    // Phase B: Construct two vertex-disjoint paths from the sink to the source.
    // ****************************************************

    int startDFS = 0;
    int maxReachedOnPath = 0;               // Highest reached vertex on first path.
    int maxReachedFrom;                     // Vertex from which the highest reached vertex was encountered.
    std::vector<int> vertsForDFS;           // Stack of vertices to explore for depth first search
    
    while (true) {
        vertsForDFS.push_back(startDFS);
        verts[startDFS].ReachedInPhaseB();

        // Phase B.1.  Try to reach higher up the existing first path by a disjoint path.

        while (!vertsForDFS.empty()) {
            int exploreVert = vertsForDFS.back();       // Next vertex for DFS
            vertsForDFS.pop_back();
            // Loop over the predecessors of exploreVert
            int predIdx = predIndex[exploreVert];
            int numPreds = predIndex[exploreVert + 1] - predIdx;
            for (int i = 0; i < numPreds; i++, predIdx++) {
                int predVert = predecessors[predIdx];             // Next predecessor of exploreVert
                if (predVert <= i && predVert >= N) {
                    assert("Invalid predecessor" && false);
                    return -2;
                }
                if (verts[predVert].TestOnFirst()) {
                    if (predVert > maxReachedOnPath) {
                        // Just reached a new furthest vertex on existing "first" path.
                        maxReachedOnPath = predVert;
                        maxReachedFrom = exploreVert;
                    }
                    continue;               // Don't add vertex on existing first path to the DFS stack.
                }
                // Add predVert to the DFS stack if it has not already been reached.
                if (!verts[predVert].ReachedInPhaseB()) {
                    vertsForDFS.push_back(predVert);
                    verts[predVert].succ = exploreVert;
                }
            }
        }
        if (maxReachedOnPath == N - 1) {
            //
            break;          // Reached the source node -- have completed the second path.
        }

        // Phase B.2. Backtrack though DFS search path to put vertices on "first" path

        int backtrackVert = maxReachedFrom;
        while (!verts[backtrackVert].TestOnFirst()) {
            verts[backtrackVert].SetOnFirst();
            backtrackVert = verts[backtrackVert].succ;
        }

        // Phase B.3. Backtrack old "first" path, remove vertices from that path, and start new DFS.

        backtrackVert = maxReachedOnPath;
        while (true) {
            int nextBacktrack = verts[backtrackVert].succ;
            if (verts[nextBacktrack].TestReachedPhaseB()) {
                break;
            }
            verts[nextBacktrack].ResetOnFirst();
            backtrackVert = nextBacktrack;
        }
        if (backtrackVert == maxReachedOnPath) {
            return -1;                              // Failed to find path onward. There are not two vertex-disjoint paths
        }

        verts[maxReachedOnPath].succ = maxReachedFrom;   // Link top part of path to new lower part of path
        startDFS = backtrackVert;
        maxReachedOnPath = -1;
    }

    // Phase B.4.  Label the nodes on the second (right) path as being on the second path.
    //             Currently all nodes on both paths are labelled as being on the first path.

    verts[N - 1].SetOnSecond();
    for (curVert = maxReachedFrom; curVert != 0; curVert = verts[curVert].succ) {
        verts[curVert].SetOnSecond();
        verts[curVert].ResetOnFirst();
    }
    verts[0].SetOnSecond();

    // *****************************************************
    // PHASE C.  Gather preliminary information on
    //      (a) bottleneck points (implication points) for the
    //          path 1 and path 2.
    // This identifies potential two-vertex bottlenext vertices.
    // *****************************************************

    std::vector<pivotInfo> leftPathPivots;
    std::vector<pivotInfo> rightPathPivots;

    for (int leftFlag = 0; leftFlag <= 1; leftFlag++) {
        std::vector<pivotInfo>& thisPivot = (leftFlag == 0) ? leftPathPivots : rightPathPivots;
        int maxReachSamePath = 0;
        int maxReachOther = 0;
        int startVert = 0;
        // (There is duplicate work done for the first DFS on the two paths. Is it worth avoiding?)
        while (true) {
            // Do a DFS for reachable verts. 
            // Stop whenever reach a vertex on either path.
            vertsForDFS.push_back(startVert);
            verts[startVert].ReachedInPhaseC(leftFlag);
            while (!vertsForDFS.empty()) {
                int exploreVert = vertsForDFS.back();       // Next vertex for DFS
                vertsForDFS.pop_back();
                // Loop over the predecessors of exploreVert
                int predIdx = predIndex[exploreVert];
                int numPreds = predIndex[exploreVert + 1] - predIdx;
                for (int i = 0; i < numPreds; i++, predIdx++) {
                    int predVert = predecessors[predIdx];             // Next predecessor of exploreVert
                    int onEither = verts[predVert].TestOnPath();
                    if (!onEither) {
                        vertsForDFS.push_back(predVert);
                    }
                    else {
                        // Update the appropriate max reached values
                        if (verts[predVert].TestOnPath(leftFlag)) {
                            maxReachSamePath = std::max(maxReachSamePath, predVert);
                        }
                        if (verts[predVert].TestOnPath(1 - leftFlag)) {
                            maxReachOther = std::max(maxReachOther, predVert);
                        }
                        if (predVert != N - 1) {
                            // Walk down the path, adding new verts to the DFS.
                            while (true) {
                                predVert = verts[predVert].succ;
                                if (verts[predVert].ReachedInPhaseC(leftFlag)) {
                                    break;
                                }
                                vertsForDFS.push_back(predVert);
                            }
                        }
                        else {
                            vertsForDFS.clear();        // Reached N-1 (source vertex). Stop the DFS.
                            break;                      // Break out of for loop, to immediately stop the DFS.
                        }
                    }
                }
            }
            thisPivot.push_back(pivotInfo(startVert, maxReachOther));
            startVert = maxReachSamePath;
            if (startVert == N - 1) {
                break;
            }
        }
        thisPivot.push_back(pivotInfo(N - 1, N - 1));
    }

#define SkipD false
    if (SkipD) {
        return -4;
    }

    // *****************************************************
    // PHASE D. Create the lists of pair groups.
    // *****************************************************

    // The PairGroupsInfo array data uses negative entries to hold
    //   information about a group as it is being formed:
    // PairGroupInfo[4*g]   equals -1 if the list of vertices on the left side is still being formed.
    // PairGroupInfo[4*g+1] equals -1 if the list of vertices on the right side is still being formed.
    // PairGroupInfo[4*g+2] equals -K if the list of the vertice on the left side has no members yet, 
    //                                but will start get members once vertex K is reached in left path.
    // PairGroupInfo[4*g+3] equals -K if the list of the vertice on the right side has no members yet, 
    //                                but will start get members once vertex K is reached in right path.
    // Postive values for these have the same meaning as in the arrays when returned.
    // It can happen that a group gets left members but not right members (or vice-versa). In this
    //     case, the group will be removed from the list.
    
    // Current number of groups, either completed and undergoing formation.
    // One side may have multiple active groups. The other side can have at most one active group.
    // "Active" means that the group may still get more members on that side.
    int gpCount = 0;
    int firstPendingGroupLeft = INT_MAX;
    int firstPendingGroupRight = INT_MAX;

    // Walk up the pivot paths, creating pair groups.
    // The vertices for the pairs come directly from the leftPathPivots and rightPathPivots.
    // Lists of Left Pivot Points and Right Pivot Points are traversed in vertex order.
    int leftPivotIndex = 0;     // Pointer to the next entry in the leftPathPivots
    int rightPivotIndex = 0;    // Pointer to the next entry in the rightPathPivots
    int prevReachedRight = 0;
    int prevReachedLeft = 0;
    while (true) {
        assert(gpCount == NumPairGroups());
        int nextLeftPivot = leftPathPivots[leftPivotIndex].PivotVert;
        int nextRightPivot = rightPathPivots[rightPivotIndex].PivotVert;
        int leftFirst = (nextLeftPivot < nextRightPivot) ? 1 : 0;

        // Handle left and right dually.  (Choose values and *references* appropriately for the duality.)
        int& thisSidePivot = leftFirst ? nextLeftPivot : nextRightPivot;
        if (thisSidePivot == N - 1) {
            PairGroupsInfo.resize(PairGroupsInfo.size() - 4); // Discard the final group (which involves just N-1)
            break;
        }
        std::vector<int>& GroupMembersThisSide = leftFirst ? GroupMembersLeft : GroupMembersRight;
        int reachedOtherPath = leftFirst ? leftPathPivots[leftPivotIndex].MaxReachOtherPath : rightPathPivots[rightPivotIndex].MaxReachOtherPath;
        int& prevReachedOther = leftFirst ? prevReachedLeft : prevReachedRight;
        int& firstPendingThisSide = leftFirst ? firstPendingGroupLeft : firstPendingGroupRight;
        int& firstPendingOtherSide = (1 - leftFirst) ? firstPendingGroupLeft : firstPendingGroupRight;
        //assert(thisSidePivot != N - 1);
        int& thisSidePivotIndex = leftFirst ? leftPivotIndex : rightPivotIndex;
        if (thisSidePivot > 0) {
            GroupMembersThisSide.push_back(thisSidePivot);      // Add the pivot index to the list of  DIP vertices
        }
        thisSidePivotIndex++;
        assert(thisSidePivotIndex - 1 == GroupMembersThisSide.size());

        // BUGGY AS HELL ....

        if (prevReachedOther < reachedOtherPath) {
            // Have an essential new crossing from this side to the other side.
            prevReachedOther = reachedOtherPath;
            // For all open groups, mark this side as complete
            // Except in the case of crossing edges, where the top group starts at the next vertex on this side.
            int g = gpCount - 1;
            bool crossingEdges = false;
            while (g >= 0 && PairGroupsInfo[4 * g + (1 - leftFirst)] == -1) {
                int& startGp = PairGroupsInfo[4 * g + 2 + (1 - leftFirst)];
                assert(PairGroupsInfo[4 * g + 2 + (1 - leftFirst)] >= 0 || g == gpCount - 1);
                if (startGp >= 0) {
                    int groupSizeThisSide = GroupMembersThisSide.size() - startGp;
                    assert(groupSizeThisSide > 0);
                    PairGroupsInfo[4 * g + (1 - leftFirst)] = groupSizeThisSide;  // Completes this side of the group.
                }
                else if ( thisSidePivot < -startGp ) {
                    // Found a the second edge of a pair of crossing edges.
                    // Start adding members to the group at this one. 
                    assert(g == gpCount - 1);
                    assert(-PairGroupsInfo[4 * (gpCount - 1) + 2 + (1 - leftFirst)] ==
                        (leftFirst ? leftPathPivots[leftPivotIndex].PivotVert : rightPathPivots[rightPivotIndex].PivotVert));  // No gap
                    assert(leftFirst ? (rightPivotIndex - 1 == PairGroupsInfo[4 * (gpCount - 1) + 3])
                        : (leftPivotIndex - 1 == PairGroupsInfo[4 * (gpCount - 1) + 2]));
                    assert(((leftFirst ? leftPathPivots[leftPivotIndex].PivotVert : rightPathPivots[rightPivotIndex].PivotVert))
                        == -PairGroupsInfo[4 * (gpCount - 1) + 2 + (1 - leftFirst)]);
                    assert (firstPendingThisSide == gpCount - 1);
                    crossingEdges = true;
                    PairGroupsInfo[4 * (gpCount - 1) + 2 + (1 - leftFirst)] = GroupMembersThisSide.size();
                }
                else {
                    // This is a vertex with incoming and outgoing edges. 
                    // Set the top group to have size one on this side.
                    PairGroupsInfo[4 * (gpCount - 1) + 2 + (1 - leftFirst)] = GroupMembersThisSide.size() - 1;
                    PairGroupsInfo[4 * (gpCount - 1) + (1 - leftFirst)] = 1;
                }
                g--;
            }
            firstPendingThisSide = INT_MAX;
            if (gpCount==0 || !crossingEdges) {
                // This is not the second of a pair of crossing edges, nor one of the starting (zero) pivots.  
                // Add a new group.
                PairGroupsInfo.resize(4 * gpCount + 4);
                PairGroupsInfo[4 * gpCount] = PairGroupsInfo[4 * gpCount + 1] = -1;
                PairGroupsInfo[4 * gpCount + 2 + (1 - leftFirst)] = GroupMembersThisSide.size();
                PairGroupsInfo[4 * gpCount + 2 + leftFirst] = -reachedOtherPath;
                // Last one is negative since the vertex is not reached on other path yet.
                firstPendingOtherSide = std::min(gpCount, firstPendingOtherSide);
                gpCount++;
            }
        }
        // Open this side as the second side of a group if appropriate.
        assert( firstPendingThisSide==INT_MAX ||
            (thisSidePivot >= -PairGroupsInfo[4 * firstPendingThisSide + 2 + (1 - leftFirst)])
            == (thisSidePivot == -PairGroupsInfo[4 * firstPendingThisSide + 2 + (1 - leftFirst)]));
        if (firstPendingThisSide < gpCount && thisSidePivot!=0 &&
                    thisSidePivot == -PairGroupsInfo[4 * firstPendingThisSide + 2 + (1 - leftFirst)]) {
            PairGroupsInfo[4 * firstPendingThisSide + 2 + (1 - leftFirst)] = GroupMembersThisSide.size();
            firstPendingThisSide = (firstPendingThisSide == gpCount - 1) ? INT_MAX : firstPendingThisSide;
            assert(firstPendingThisSide == INT_MAX || PairGroupsInfo[4 * firstPendingThisSide + 2 + (1 - leftFirst)] < 0);
        }

    }

    return NumPairGroups();
 }