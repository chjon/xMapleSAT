
//
// TwoVertexBottlenecks.cpp
//   Version 1.0, January 13, 2023
//   Version 1.1. October 17, 2023. Minor revisions & Phase C performance fix.
//   Version 2.0. October 30, 2023. Major revision, bug fix, and interface change.
// .
//   Author: Sam Buss (sbuss@ucsd.edu) with Jonathan Chung, Vijay Ganesh and Albert Oliveras
//   This code has no warranties of correctness or appropriateness.
//   May be used freely, however acknowledgement of use 
//         is expected and appreciated.
// 
//  ************************************************
//  * Main routines for TwoVertexBottlenecks       *
//  *                                              *
//  * See TwoVertexBottleneck.h for usage.         *
//  ************************************************

#include "TwoVertexBottlenecks.h"
#include <algorithm>
#include <limits>

inline void UpdateMax(int& runningMax, int newValue) {
    if (runningMax < newValue) {
        runningMax = newValue;
    }
}

inline void UpdateMin(int& runningMin, int newValue) {
    if (runningMin > newValue) {
        runningMin = newValue;
    }
}

// ********************************************
// Internal data structures.
// ********************************************

//   vertInfo-- Used in Phases A, B, C
//   Tracks for each vertex: whether it is on a path
//           and whether it has been visited the depth first searches.

class vertInfo {
public:
    int status = 0;     // Status
    int succ = -1;      // Successor node, i.e, on the path towards the sink (vertex 0).

    // Handle membership on first path.  
    bool TestOnFirst() const { return status & 0x01; }  // Test whether on first path.
    void SetOnFirst() { status = status | 0x01; }       // Mark as on first path
    void ResetOnFirst() { status = (status & ~0x01); }  // Mark as not on first path

    bool TestOnSecond() const { return status & 0x02; } // Test whether on second path
    void SetOnSecond() { status = status | 0x02; }      // Mark as on second path

    // TestOnPath returns 0, 1, 2, 3 - If on neither path, path one, path two or both.
    int TestOnPath() const { return (status & 0x03); }   
    // Test if on path two (leftFlag==1) or on path 1 (leftFlag==0)
    bool TestOnPath(int leftFlag) const { return status & (0x01 << leftFlag); } 

    bool TestReachedPhaseB() const { return status & 0x04; }
    // Mark as reached in Phase B.  Returns true if already was reached in Phase B
    bool ReachedInPhaseB() {
        bool ret = status & 0x04;
        status = status | 0x04;
        return ret;
    }
};

class reachInfo {
    // maxDirectReachOnPath - max. reachable node on the path.
    // minAncestorOnPath - min reachable node on the path,
    //                       allowing traversing vertices on the *other* path.
public:
    reachInfo(int maxReach, int minAncestor)
        : maxDirectReachOnPath(maxReach), minAncestorOnPath(minAncestor) {}
    int maxDirectReachOnPath;
    int minAncestorOnPath;
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

class pivotInfoPlus {
    // In the arrays of left/rightPathPivot's:
    //     First PivotVert 0. (Sink)   Last PivotVert is N-1. (Source)
    //     MaxReachOtherPath may not be non-decreasing, but it
    //         will be once points on the other path which are not pivots are skipped.
    //     MinAncestorOtherPath is non-decreasing.
    //     MaxReachOtherPath is initially used for the max reachable vertex on the same path,
    //         but then is replaced with the max reachable vertex on the other path.

public:
    pivotInfoPlus() { assert(false); }
    pivotInfoPlus(int vert, int maxReach) 
    : PivotVert(vert), MaxReachOtherPath(maxReach) {}

    int PivotVert;
    int MaxReachOtherPath;
    int MinAncestorOtherPath;
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
            // Reached the source node -- have completed the second path.
            break;
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
            SingleVertBottleneck = maxReachedOnPath;    // First bottlenext vertex.
            return -1;                              // Failed to find path onward. There are not two vertex-disjoint paths
        }

        verts[maxReachedOnPath].succ = maxReachedFrom;   // Link top part of path to new lower part of path
        startDFS = backtrackVert;
        maxReachedOnPath = -1;
    }

    // Phase B.4.  Label the nodes on the second (right) path as being on the second path.
    //             Currently all nodes on both paths are labelled as being on the first path.

    // The vertex (after N-1) at the top of path 1 is verts[N-1].succ.
    // The vertex (after N-1) at the top of path 2 is lastVertSecondPath.
    // After that, the verts[].succ pointers give successors along the paths.
    const int& lastVertSecondPath = maxReachedFrom;

    // The parts of the code below setting infomation abut PathA and PathB are
    //     for reporting purposes only. The calling program might find
    //     it useful to have information about PathA and PathB, but this 
    //     is not needed anymore by CalcBottlenecks.
    PathA.clear();
    verts[N - 1].SetOnSecond();
    PathA.push_back(N - 1);
    for (curVert = lastVertSecondPath; curVert != 0; curVert = verts[curVert].succ) {
        verts[curVert].SetOnSecond();
        verts[curVert].ResetOnFirst();
        PathA.push_back(curVert);
    }
    verts[0].SetOnSecond();
    PathA.push_back(0);

    PathB.clear();
    for (curVert = N - 1; ; curVert = verts[curVert].succ) {
        assert(verts[curVert].TestOnFirst());
        PathB.push_back(curVert);
        if (curVert == 0) {
            break;
        }
    }


    // *****************************************************
    // PHASE C-1.  Gather preliminary information on
    //       reachability from path nodes to other path nodes.
    //   Paths 1 and 2 (right path and left path, respectively) are
    //   now fixed and will not change.
    //   We gather information for each vertex x on path 1 (resp. path 2)
    //     (a) the maximum nodes on paths 1 and 2 which
    //         can be reached from x without first encountering any
    //         vertex from *either* path ("maxDirectReachOnPath")
    //     (b) the minimum node on path 2 (resp. path 1) that is an
    //         predecessor of x (that is, on a path towards N-1).
    //         This node ("minAncestorOnPath") is on the *other* path from x's path.
    //   To do this, we collect these data items for *every* vertex
    //     in the graph, obtaining the values for min-reachable and
    //     max-directly-reachable nodes in both paths.
    //   These data items are saved in DirectReach1 and DirectReach2,
    //     in reverse order (i.e., from N-1 down to 0).
    //   DirectReach1 and DirectReach2 for reachability *to* vertices on paths 1 and 2.
    //   We are interested in these values for all vertices in *both* paths.
    // *****************************************************

    std::vector<reachInfo> directReach1;
    std::vector<reachInfo> directReach2;

    directReach1.emplace_back(N - 1, N - 1);  // First entry is for vertex N-1 (the source)
    directReach2.emplace_back(N - 1, N - 1);  // ditto
    int lastPredIdx = predIndex[N - 1] - 1;  // No predecessor is traversed yet.
    for (int i = N - 2; i >= 0; i--) {
        // Traverse sequentially backwards (from source to sink)
        //   through all vertices in the graph.
        int firstPredIdx = predIndex[i];        // Index in predecessors[] of the first precessor of vertex i.
        int newMin1 = N - 1, newMin2 = N - 1;
        int newMax1 = 0, newMax2 = 0;
        while (lastPredIdx >= firstPredIdx) {
            // lastPredIdx is an index for a predecessor of thisVert. Handle this now.
            // Decrement lastPredIdx for loop control.
            int predVertNum = predecessors[lastPredIdx--];
            int testOnPaths = verts[predVertNum].TestOnPath();
            assert(testOnPaths != 3 || predVertNum == N - 1);
            if (testOnPaths == 0) {     // If on neither path
                UpdateMax(newMax1, directReach1[N - 1 - predVertNum].maxDirectReachOnPath);
                UpdateMax(newMax2, directReach2[N - 1 - predVertNum].maxDirectReachOnPath);
            }
            else {
                if (testOnPaths & 0x01) {     // If predecessor is on path one
                    UpdateMax(newMax1, predVertNum);
                }
                if (testOnPaths & 0x02) {     // If predecessor is on path two (only on both if vertex N-1).
                    UpdateMax(newMax2, predVertNum);
                }
            }
            UpdateMin(newMin1,
                testOnPaths == 1 ? predVertNum : directReach1[N - 1 - predVertNum].minAncestorOnPath);
            UpdateMin(newMin2,
                testOnPaths == 2 ? predVertNum : directReach2[N - 1 - predVertNum].minAncestorOnPath);
        }
        directReach1.emplace_back(newMax1, newMin1);
        directReach2.emplace_back(newMax2, newMin2);
    }

    assert((directReach1.back().maxDirectReachOnPath == N - 1) == (directReach2.back().maxDirectReachOnPath == N - 1));
    if (directReach1.back().maxDirectReachOnPath == N - 1) {
        return 0;           // Three connected. Hence no two vertex bottlenecks!
    }

    // ++++++++++++++++++++++++++++++++++++++++++++++
    // Phase C-2.  
    // Build the lists of vertices for paths 1 and 2
    //   that could be in a two-vertex bottleneck,
    //   i.e., removing vertices that are bypassed.
    // At the same time, preserve the information about
    //    the min-reachable on the other path,
    // And update the max-direct-reachable vertex on the other
    //    path to include *only* reachable from nodes which are
    //    *below* the current node. 
    //    (This makes the max-reachable-other values non-decreasing)
    // All this information is stored in leftPathPivots
    //    and rightPathPivots
    // **********************************************
    std::vector<pivotInfoPlus> leftPathPivots;
    std::vector<pivotInfoPlus> rightPathPivots;
    for (int leftFlag = 0; leftFlag <= 1; leftFlag++) {
        // leftFlag==0 - finding the pivots for the left path
        // leftFlag==1 - finding the pivots for the right path.
        std::vector<pivotInfoPlus>& thisPathPivots = (leftFlag == 0) ? leftPathPivots : rightPathPivots;
        const std::vector<reachInfo>& directReachThis = (leftFlag == 0) ? directReach2 : directReach1;
        const std::vector<reachInfo>& directReachOther = (leftFlag == 0) ? directReach1 : directReach2;
        int nextOnPath = leftFlag == 0 ? lastVertSecondPath : verts[N - 1].succ;
        do {
            thisPathPivots.emplace_back(nextOnPath, directReachThis[N - 1 - nextOnPath].maxDirectReachOnPath);
            nextOnPath = verts[nextOnPath].succ;
        } while (nextOnPath != 0);
        // Reverse the elements in the list of potential pivots
        // This now holds the entire current path, except not vertex 0.
        std::reverse(thisPathPivots.begin(), thisPathPivots.end());
        // Scan the list of potential pivots and remove the ones that are jumped past
        //    Removal is governed by maxDirectReachOnPath value pointing higher in the path.
        // In addition set the MinAncestorOtherPath and MaxReachOtherPath values
        //    (The MaxReachOtherPath data replaces the maxDirectReachOnPath value.)
        auto toEntry = thisPathPivots.begin();
        auto fromEntry = toEntry;  
        int runningMaxAncestorOnPath = directReachThis[N - 1].maxDirectReachOnPath;  // Pick up sink's (Vertex 0's) direct reach
        int runningMaxReachOtherPath = directReachOther[N - 1].maxDirectReachOnPath; // Ditto
        while (fromEntry != thisPathPivots.end()) {
            int thisPivotVert = (fromEntry++)->PivotVert;  // Increment fromEntry for loop control.
            if (thisPivotVert >= runningMaxAncestorOnPath) {
                toEntry->PivotVert = thisPivotVert;
                toEntry->MaxReachOtherPath = runningMaxReachOtherPath;
                (toEntry++)->MinAncestorOtherPath = directReachOther[N - 1 - thisPivotVert].minAncestorOnPath;
            }
            UpdateMax(runningMaxReachOtherPath, directReachOther[N - 1 - thisPivotVert].maxDirectReachOnPath);
            UpdateMax(runningMaxAncestorOnPath, directReachThis[N - 1 - thisPivotVert].maxDirectReachOnPath);
        }
        thisPathPivots.resize(std::distance(thisPathPivots.begin(), toEntry));
        if (thisPathPivots.size() == 0) {
            return 0;                   // No potential DIP members. The graph is three connected.
        }
        thisPathPivots.emplace_back(N - 1, N - 1);
    }

    // ****************************************************
    // Phases D-1 & D-2 - Complete the analysis of two-vertex bottlenecks.
    //   Fill in the information in VertListA and VertListB with
    //   the complete specification of allowable two vertex bottlenecks
    //   the minAncestor information.
    // 
    // Phase D-1.  Compute for each vertex x in each path (left and right), 
    //   an upper and lower bound on the vertices y in the other path (right and left, resp.)
    //   such that (x,y) forms a two vertex bottleneck (TVB) pair.
    //   These upper and lower bounds may not actually be present as vertices 
    //   in the other path.
    //   At the same, remove all vertices that do not belong to a valid TVB pair.
    // 
    //   The lower and bound values are stored in the fields minPairIdx and maxPairIdx.
    //   Phase D-2 will overwrite these values with indices from the other vertList
    //   array.
    // ****************************************************
    for (int leftFlag = 0; leftFlag <= 1; leftFlag++) {
        // leftFlag==0 - finding the information for the left path (VertListA)
        // leftFlag==1 - finding the information for the right path (VertListB)
        const std::vector<pivotInfoPlus>& thisPathPivots = (leftFlag == 0) ? leftPathPivots : rightPathPivots;
        const std::vector<pivotInfoPlus>& otherPathPivots = (leftFlag == 0) ? rightPathPivots : leftPathPivots;
        std::vector<VertPairInfo>& vertListX = (leftFlag == 0) ? VertListA : VertListB;
        int butLastOtherPathPivotVert = (otherPathPivots.end() - 2)->PivotVert;
        auto otherPathPivot = otherPathPivots.begin();
        auto curPivot = thisPathPivots.begin();
        for ( ; curPivot->PivotVert != N - 1; curPivot++) {
            while (otherPathPivot->MaxReachOtherPath <= curPivot->PivotVert) {
                otherPathPivot++;
            }
            if (otherPathPivot == otherPathPivots.begin()) {
                continue;    // Vertex 0 reached beyond curPivot
            }
            int upperBdPair = (otherPathPivot-1)->PivotVert;
            int lowerBdPair = curPivot->MaxReachOtherPath;   // Max reachable by a vertex below the curPivot vertex.
            if (lowerBdPair <= upperBdPair && (upperBdPair < N - 1 || lowerBdPair <= butLastOtherPathPivotVert)) {
                vertListX.emplace_back(curPivot->PivotVert, lowerBdPair, upperBdPair, curPivot->MinAncestorOtherPath);
            }
        }
    }
    assert(VertListA.empty() == VertListB.empty());
    if (VertListA.empty()) {
        return 0;    // No two vertex bottlenecks - Source and sink are three connected
    }

    // ****************************************************
    // Phase D-2. Replace the upper and lower bounds from D-1 on
    //     vertices that form TVB pairs with indices into the other vertList.
    //     The result is that vertex VertListA[i] and VertListB[j] form
    //     a two vertex bottleneck (TVB) pair iff
    //       VertListA[i].minPairIdx <= j <= VertListA[i].maxPairIdx.
    //     and iff
    //       VertListB[j].minPairIdx <= i <= VertListB[j].maxPairIdx.
    // ****************************************************
    for (int leftFlag = 0; leftFlag <= 1; leftFlag++) {
        // leftFlag==0 - computing the information for the left path (VertListA)
        // leftFlag==1 - computing the information for the right path (VertListB)
        std::vector<VertPairInfo>& vertListThis = (leftFlag == 0) ? VertListA : VertListB;
        std::vector<VertPairInfo>& vertListOther = (leftFlag == 0) ? VertListB : VertListA;
        int minPairI = 0;
        int maxPairI = 0;
        int s = (int)vertListOther.size();
        for (VertPairInfo& vert : vertListThis) {
            int minpr = vert.minPairIdx;   // The min **vertex number** that can be paired
            while (vertListOther[minPairI].vertNum < minpr) {
                minPairI++;
            }
            vert.minPairIdx = minPairI; // Set the **index** of the minimum vertex to pair
            int maxpr = vert.maxPairIdx;   // The max **vertex number** that can be paired
            while (maxPairI < s - 1 && vertListOther[maxPairI + 1].vertNum <= maxpr) {
                maxPairI++;
            }
            vert.maxPairIdx = maxPairI; // Set the **index** of the maximum vertex to pair
        }
    }
    return 1;       // DIP's were found.

 }