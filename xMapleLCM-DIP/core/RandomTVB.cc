//
// RandomTVB.cpp
//   Version 1.0, Started 11/24/2024
// .
//   Author: Sam Buss (sbuss@ucsd.edu) with Jonathan Chung, Vijay Ganesh and Albert Oliveras
//   This code has no warranties of correctness or appropriateness.
//   May be used freely, however acknowledgement of use 
//         is expected and appreciated.
// 

#include "RandomTVB.h"
#include "TwoVertexBottlenecks.h"
#include <vector>


bool ChooseRandomTVB(
    const std::vector<TwoVertexBottlenecks::VertPairInfo>& vertListA,
    const std::vector<TwoVertexBottlenecks::VertPairInfo>& vertListB,
    int& vertIdxA, int& vertIdxB, bool avoidInit) {

    // Get the total number of DIPs
    int numDIPs = 0;
    for (const TwoVertexBottlenecks::VertPairInfo& vpiA : vertListA) {
        numDIPs += vpiA.maxPairIdx - vpiA.minPairIdx + 1;
    }

    int avoidNum = (int)avoidInit;
    if (numDIPs <= avoidNum) {
        return false;
    }

    int randNum = RandGenTVB(numDIPs - avoidNum);
    assert(randNum >= 0 && randNum < numDIPs - avoidNum);
    int randomDIP = randNum + avoidNum;

    for (int i = 0; ; i++) {
        const TwoVertexBottlenecks::VertPairInfo& vpiA = vertListA[i];
        int numPairs = vpiA.maxPairIdx - vpiA.minPairIdx + 1;
        if (numPairs > randomDIP) {
            vertIdxA = i;
            vertIdxB = vpiA.minPairIdx + randomDIP;
            return true;
        }
        randomDIP -= numPairs;
    }
}

