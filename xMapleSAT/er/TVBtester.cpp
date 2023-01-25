
// *******************************
// Tester coder for TwoVertexBottleNecks.
// *******************************

#include "TwoVertexBottlenecks.h"

void SimpleTests();

int main() {
    
    SimpleTests();

}

void SimpleTests()
{
    TwoVertexBottlenecks myTVB;
    int i;

    // Test "A": Simple path 0 - 2 - 4 - 6.  No second path.
    int N_A = 7;
    int preds_A[] =   { 2,3,4,5,6,6 };
    int predIdx_A[] = { 0,1,2,3,4,5,6 };
    i = myTVB.CalcBottlenecks(N_A, preds_A, predIdx_A);
    assert(i == -1);    // No second vertex-disjoint path

    // Test "B":  Like Test "A", but add edges 0 - 1 and 3 - 4 
    int N_B = 7;
    int preds_B[] =   { 2,1, 3, 4, 4,5, 6, 6 };
    int predIdx_B[] = { 0,   2, 3, 4,   6, 7, 8 };
    i = myTVB.CalcBottlenecks(N_B, preds_B, predIdx_B);
    assert(i == 2);    

    // Test "C":  Like Test "B", but with edge  4-5 instead of 3-5. Makes path reroute in Phase B. 
    int N_C = 7;
    int preds_C[] =   { 2,1, 3, 4, 4,   6,5,  6 };
    int predIdx_C[] = { 0,   2, 3, 4,   5,    7,  8 };
    i = myTVB.CalcBottlenecks(N_C, preds_C, predIdx_C);
    assert(i == -1);    // No second vertex-disjoint path

    // Test "D":  Like Test "D", but adding edge 2-5. Finds 2 paths, after re-routing. 
    int N_D = 7;
    int preds_D[] =   { 2,1, 3, 4,5, 4,   6,5,  6 };
    int predIdx_D[] = { 0,   2, 3,   5,   6,    8,  9 };
    i = myTVB.CalcBottlenecks(N_D, preds_D, predIdx_D);
    assert(i == 2);

    // Test "E":  Like Test "D", but adding edges 4-7,6-7,6-8,7-8. Finds 2 paths, after three times re-routing. 
    int N_E = 9;
    int preds_E[] = { 2,1, 3, 4,5, 4,   6,5,7,  6,  7,8,  8 };
    int predIdx_E[] = { 0,   2, 3,   5,   6,      9,  10,   12, 13 };
    i = myTVB.CalcBottlenecks(N_E, preds_E, predIdx_E);
    assert(i == 3);

    // Test "F":  Like Test "E", but switches out-order of node 6. Finds 2 paths, after twice re-routing. 
    int N_F = 9;
    int preds_F[] =   { 2,1, 3, 4,5, 4,   6,5,7,  6,  8,7,  8 };
    int predIdx_F[] = { 0,   2, 3,   5,   6,      9,  10,   12, 13 };
    i = myTVB.CalcBottlenecks(N_F, preds_F, predIdx_F);
    assert(i == 3);

    // Test "G"
    int N_G = 12;
    int preds_G[] =   { 2,1, 3,4, 4,6, 5, 5, 10,9,6, 7, 8,9, 10, 11, 11 };
    int predIdx_G[] = { 0,   2,   4,   6, 7, 8,      11,12,  14, 15, 16, 17 };
    i = myTVB.CalcBottlenecks(N_G, preds_G, predIdx_G);
    assert(i == 3);
}
