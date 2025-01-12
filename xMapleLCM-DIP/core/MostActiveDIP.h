#pragma once
// *********************************************************
// MostActiveDip.h -- 
//    Version 1.0 - started November 2, 2024 
// 
//     Author: Sam Buss (sbuss@ucsd.edu) with Jonathan Chung, Vijay Ganesh and Albert Oliveras
//     This code has no warranties of correctness or appropriateness.
//     May be used freely, however acknowledgement of use
//         is expected and appreciated.
//
//  These routines take as input a "VertListA" and "VertListB" as output
//     by the routine CalcBottlenexts(...) in TwoVertexBottlenecks.cpp,
//  plus a function computing activity values for the vertices in VertListA and
//     and VertListB
//  and alsp a function combining two activity values to produce a composite
//     activity level, e.g., by summing or multiplying activity levels.
// 
//  It returns the two-vertex bottleneck consisting of nodes (a,b)
//     such that the sum/product/other combination of the activity levels 
//     of vertex a and vertex b is maximized.
// *********************************************************

#ifndef MOST_ACTIVE_DIP_H
#define MOST_ACTIVE_DIP_H

#include <vector>
#include "TwoVertexBottlenecks.h"

// *********************************************************
//  To compile these routines, you must typedef ActivityValue_f to be the 
//  correct data type for your application.  In addition, the
//  functions
// 
//     template<typename ACTIVITY>
//     ACTIVITY CombineActivities(ACTIVITY aValue, ACTIVITY bValue);
// and
//     template<typename ACTIVITY>
//     ACTIVITY ActivityOf(int vertNum);
// 
//  must be defined before linkage. See examples in TVBtester.cpp
//  -----------------------

// This function (user-defined somewhere else)
//   - must take as input a vertex number of a literal (namely, a vert in VertListA or VertListB)
//     and return the activity level for the literal.
template<typename ACTIVITY>
ACTIVITY ActivityOf(int vertNum);

// This function combines two ACTIVITY's into a composite ACTIVITY.
//   Usually, by adding or multiplying (whilst avoiding overflow).
template<typename ACTIVITY>
ACTIVITY CombineActivities(ACTIVITY aValue, ACTIVITY bValue);
// *********************************************************

//
// FindMaxActivityDIP & FindMaxActivityDIP_idx
// 
//   - Input is the VertListA and VertListB from a TwoVertexBottlenecks object;
//   - Return values, return_a and return_b are
//        a pair of vertices (a,b), from VertListA and VertListB, respectively
//        such that ActivityOf(a) + GetActivityOf(b) is maximized.
//     OR return_a_idx and return_b_idx are the indices in VertListA and VertListB
//        for these vertices.
//     That is, return_a is equal to VertListA[return_a_idx].vertNum (similarly for "b")
//   - The function ActivityOf() is needed for obtaining activity values.
//   - The function CombineActivities() is needed for combining pairs of activity values.
//
// These FindMaxActivityDIP... routines should be called 
//     only if the TwoVertexBottlenecks function found DIP's.

template<typename ACTIVITY>
ACTIVITY FindMaxActivityDIP_idx(
    const std::vector<TwoVertexBottlenecks::VertPairInfo>& vertListA,
    const std::vector<TwoVertexBottlenecks::VertPairInfo>& vertListB,
    int& return_a_idx, int& return_b_idx);

template<typename ACTIVITY>
ACTIVITY FindMaxActivityDIP(
    const std::vector<TwoVertexBottlenecks::VertPairInfo>& vertListA,
    const std::vector<TwoVertexBottlenecks::VertPairInfo>& vertListB,
    int& return_a, int& return_b)
{
    int return_a_idx;
    int return_b_idx;
    ACTIVITY mostActivity = FindMaxActivityDIP_idx<ACTIVITY>(vertListA, vertListB, return_a_idx, return_b_idx);
    return_a = vertListA[return_a_idx].vertNum;
    return_b = vertListB[return_b_idx].vertNum;
    return mostActivity;
}

template<typename ACTIVITY>
ACTIVITY FindMaxActivityDIP_idx(
    const std::vector<TwoVertexBottlenecks::VertPairInfo>& vertListA,
    const std::vector<TwoVertexBottlenecks::VertPairInfo>& vertListB,
    int& return_a_idx, int& return_b_idx)
{
    int sizeA = (int)vertListA.size();
    int sizeB = (int)vertListB.size();

    // bestAvailable is a FIFO queue of best activity value matches.
    //   We push and pop on the right (the top) and only pop on the left.
    // bestAvailable holds decreasing activity levels from vertListA.
    // vertListA[bestAvailableIdx[i]] gives the value where bestAvailable[i]
    //   replaces earlier (larger) bestAvailable values.
    std::vector<ACTIVITY> bestAvailable;
    std::vector<int> bestAvailableIdx;
    int curBestIdx = 0;         // Index of the leftmost entry in the bestAvailable array.
    bestAvailableIdx.push_back(0);      // bestAvailableIdx[0] := 0;
    bestAvailable.push_back(ActivityOf<ACTIVITY>(vertListA[0].vertNum));  // bestAvailable[0] activity value

    int ret_a = -1;
    int ret_b = -1;
    ACTIVITY maxValueSoFar;

    // i - index into vertListA --- also equal to one plus the largest (rightmost) bestAvailableIdx[] value.
    // j - index into vertListB
    int i = 1;
    for (int j = 0; j < sizeB; j++) {
        // Update the "bestAvailable" values:
        // First add values on the right end of the queues, overwriting (deleting) smaller values
        for (; i < sizeA; i++) {
            if (vertListB[j].maxPairIdx < i) {
                break;
            }
            // Add this value to the right end, after popping off smaller activity values.
            ACTIVITY iActivity = ActivityOf<ACTIVITY>(vertListA[i].vertNum);
            while (curBestIdx < (int)bestAvailable.size() && bestAvailable.back() < iActivity) {
                bestAvailable.pop_back();
                bestAvailableIdx.pop_back();
            }
            if (curBestIdx == (int)bestAvailable.size()) {
                bestAvailable.clear();      // Clear back to the beginning if possible.
                bestAvailableIdx.clear();
                curBestIdx = 0;
            }
            bestAvailable.push_back(iActivity);
            bestAvailableIdx.push_back(i);
        }
        // Second, remove values from the left end that are no longer available(.
        while (bestAvailableIdx[curBestIdx] < vertListB[j].minPairIdx) {
            curBestIdx++;
        }
        assert(curBestIdx < (int)bestAvailable.size());
        ACTIVITY thisActivityValue = CombineActivities<ACTIVITY>(ActivityOf<ACTIVITY>(vertListB[j].vertNum), bestAvailable[curBestIdx]);
        if (ret_a == -1 || maxValueSoFar < thisActivityValue) {
            ret_a = bestAvailableIdx[curBestIdx];
            ret_b = j;
            maxValueSoFar = thisActivityValue;
        }
    }

    return_a_idx = ret_a;
    return_b_idx = ret_b;

    return maxValueSoFar;
}

#endif // MOST_ACTIVE_DIP_H
