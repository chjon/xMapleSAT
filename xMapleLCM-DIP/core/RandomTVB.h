#pragma once
// *********************************************************
// RandomTVB.h -- 
//    Version 1.0 - started November 24, 2024 
// 
//     Author: Sam Buss (sbuss@ucsd.edu) with Jonathan Chung, Vijay Ganesh and Albert Oliveras
//     This code has no warranties of correctness or appropriateness.
//     May be used freely, however acknowledgement of use
//         is expected and appreciated.
//
//  This routine takes as input VertListA and VertListB as produced
//      by CalcBottlenecks, and returns a randomly selected two-vertex
//      bottleneck (a,b). 
//  Optionally, the minimal pair is disallowed.
//  Randomization is on the function RandGenTVB(i), which must
//      be user supplied (function prototype is in this header file).
// *********************************************************

#ifndef RANDOM_TVB_H
#define RANDOM_TVB_H

#include <vector>
#include "TwoVertexBottlenecks.h"

// ** RandGenTVB ** 
//    User supplied function, returns random integer in the range [0,n).
//    This allows the calling program to control random number generation.
int RandGenTVB(int n);      

// *********************
// ChooseRandomTVB
// Input is:
//     - the two Vert Lists from CalcBottleNecks
//     - Optional argument avoidInit: whether the randomly chosen DIP can be the 
//                  first vertices from vertListA and vertListB.
// Returns:
//     - return code is TRUE if a randomly chosen DIP is returned.
//       return code is FALSE if there is no DIP available to return.
//     - return values vertIdxA and vertIdxB are the indices in vertListA and vertListB of a TVB.
//       I.e., vertListA[vertIdxA].vertNum and vertListB[vertIdxB].vertNum
//             are the members of the two vertex bottleneck.
// *********************
bool ChooseRandomTVB(const std::vector<TwoVertexBottlenecks::VertPairInfo>& vertListA, 
                     const std::vector<TwoVertexBottlenecks::VertPairInfo>& vertListB, 
                     int& vertIdxA, int& vertIdxB, bool avoidInit = false);

#endif   // RANDOM_DIP_H
