#ifndef RLIST_H
#define RLIST_H

#include "Data.h"

RHdrNd RinitHdrNd();
RLstNd RinitLstNd(RTreeNode tnInfo);
void RinsertLstElem(RHdrNd HdrNdLst, RTreeNode tnInfo);
void RinsertLstNd(RHdrNd HdrNdLst, RLstNd LstNdElem);
Boolean RisLstEmpty(RHdrNd HdrNdLst);
RLstNd RdeleteLstElem(RHdrNd HdrNdLst, RTreeNode tnInfo);
RLstNd RdeleteLstFirst(RHdrNd HdrNdLst);
RLstNd RdeleteNextNd(RHdrNd HdrNdLst, RLstNd LstNdElem);
RNbHdr RinitNbHdr();
Boolean RisNbLstEmpty(RNbHdr nbHdrInfo);
void RinsertNeighbors(Data dataTemp, Data dataClstElem, double dist);
void RprintNbhLst(RNbHdr nbHdrInfo, DataHdr dataList1);
void RfreeNeighborhood(Data currDataPoint);
void RDeleteAll(RHdrNd hdrndlist);

#endif
