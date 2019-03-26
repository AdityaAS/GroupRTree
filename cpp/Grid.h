#ifndef GRID_H
#define GRID_H

#include "GridRTree.h"

GHdrNd constructMainGRTree(DataHdr dataList1, BCellListHd cellsList, int gMinEntries, int gMaxEntries);
void constructAuxRTrees(DataHdr dataList1, GroupListHd groupList);
void printAuxRTrees(BCellListHd cellsList);
void verifyAuxRTrees(BCellListHd cellsList);
void freeAuxRTrees(BCell cellsList);

#endif
