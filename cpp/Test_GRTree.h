#ifndef __TEST_GRTREE_H
#define __TEST_GRTREE_H

#include "GridRTree.h"

int validateNeighborhood(RNbHdr n1,RNbHdr n2 );
GHdrNd buildAuxGRTree(BCellListHd cellsList, int gMinEntries, int gMaxEntries);
void runTestCells(BCellListHd cellsList, GHdrNd GridRTree1, RHdrNd pointsRTree, DataHdr dataList1);
void runTestPoints(BCellListHd cellsList, GHdrNd GridRTree1, RHdrNd pointsRTree, DataHdr dataList1);



#endif