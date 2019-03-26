#ifndef GRIDRTREE_H

#define GRIDRTREE_H

#include "GList.h"
#include "BCell.h"
#include "RTree.h"

void insertGroupIntoGroupList(Group newGroup, GroupListHd cellsList);
void insertPointintoGroup(Group currGroup, Data currDataPoint);
Boolean GisOverLap(Region rgnRectOne, Region rgnRectTwo);
Boolean GisOverLapTotal(Region rgnRectOne, Region rgnRectTwo);
Boolean GisOverLap2(Region rgnRectOne, Region rgnRectTwo);
Region GinitRgnRect(Dimension iBottomLeft, Dimension iTopRight);
GTreeNode GinitIntNd(Dimension iBottomLeft, Dimension iTopRight);
GTreeNode GinitExtNd(Group groupElem);
void GsetRect(GLstNd lstNd, GTreeNode tnInfo);
GLstNd GpickChild(GHdrNd ptrChildLst, GTreeNode tnInfo);
Boolean GexpansionArea(Region rgnRect, GTreeNode tnInfo, Double ptrDMinExp, Region rgnNewRect);
double Garea(Region rgnRect);
Boolean GinsertTree(GHdrNd hdrNdTree, GTreeNode tnInfo, int gMinEntries, int gMaxEntries);
GHdrNd GcreateRoot(GHdrNd hdrNdTree);
void GsplitNode(GLstNd ptrChild, int gminEntries);
void GpickSeeds(GHdrNd ptrChildLst, GLstNd *lstNdChildOne, GLstNd *lstNdChildTwo);
int GisContainsCell(Region rgnRect, Group group);
Group GfindCell(GHdrNd hdrNdTree, Region rgnRect);
void GprintTree(GHdrNd hdrNdTree);
void GprintRegion(Region region);
void GgetNeighborHood(GHdrNd hdrNdTree, Data currDataPoint);
void appendNbh(Data currDataPoint, Data temp);
BCellListHd GgetCellsInRegion(GHdrNd hdrNdTree, Region regRect, Region cellRgnRect);
BCellListHd GgetCellsInRegion2(GHdrNd hdrNdTree, Region regRect, Region cellRgnRect);
BCellListHd GgetCellsInTotalRegion2(GHdrNd hdrNdTree, Region regRect, Region cellRgnRect);
int GgetCellsInRegion3(GHdrNd hdrNdTree, Region regRect, Region cellRgnRect);
unsigned int GfindOverlappingCells(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList);
Boolean GisOverLapTotal_MN(Region rgnRectOne, Region rgnRectTwo, BCell tempCell, int myrank);
unsigned int GfindOverlappingCellsTotalOverlap_MN(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList, BCell tempCell, int myrank);
int GgetCellsInRegion_MN(GHdrNd hdrNdTree, Region regRect, Region cellRgnRect, BCell tempCell, int myrank);
unsigned int GfindOverlappingCellsTotalOverlap(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList);
unsigned int GfindOverlappingCellsOptimized(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList, Region halfEpsExtended);
BCellListHd GgetCellsInRegionOptimized(GHdrNd hdrNdTree, Region regRect, Region cellRgnRect, Region halfEpsExtended);
unsigned int GfindOverlappingCells2(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList);
unsigned int GfindTotalOverlappingCells2(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList);
unsigned int GfindOverlappingCells3(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList,BCellListHd cellsList2);
GHdrNd populateAuxGridRTree(BCellListHd cellsList, int gMinEntries, int gMaxEntries);
void findReachableGroupsofGroupG(GHdrNd GRTree, Group G);
void findReachableGroupsRTree(GHdrNd GRTree, vector<Group>& groupList);
GHdrNd populateGridRTree(DataHdr dataHdrLst, GroupListHd cellsList, int gMinEntries, int gMaxEntries);
void insertPointintoBCell(BCell currBcell, Data currDataPoint);
void insertBCellIntoBCellList(BCell newBCell,BCellListHd cellsList);
void insertPointintoCellDataList(CellDataHd currCellDataList, Data currDataPoint);
Region createCellRegOfPoint(Data currDataPoint, double eps);
Region createCellRegOfDataPoint(DataPoint currDataPoint, double eps);
Region createRegionofCell(BCell bCellElem);
GHdrNd insertGroupIntoRTree(GHdrNd hdrNdTree, Group groupNode, int gMinEntries, int gMaxEntries);
void printMinGridSize();
void printMaxGridSize();
void printCellsList(BCellListHd cellsList);
void printNoOfCoreCells(BCellListHd cellsList);
void printCell(BCell bCellElem);
void printCellDataList(CellDataHd currCellDataList);
void printCellData(CellData currCellData);
Region getEpsExtendedRegion(Region cellRgnRect, double eps);
Region getEpsOptimalExtendedRegion(Region cellRgnRect, double eps);
Region getEpsExtendedRegionPoint(Data dataPoint, double eps);
void freeGRTree(GHdrNd hdrNdTree);
GHdrNd populateCellGridRTree(CellDataHd cellDataList, BCellListHd cellsList, BCell thisBCell, int gMinEntries, int gMaxEntries);
void isCorrectGRTree(GHdrNd hdrNdTree);
vector<int> findNeighbours(DataHdr dataHdrLst, int point_id, double eps);

//RHdrNd getAuxRTreeofNeighbors(BCell bCellElem, GHdrNd hdrNdTree);

#endif