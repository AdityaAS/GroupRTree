#include "BCell.h"
#include "RList.h"
#include "GridRTree.h"
#include "Grid.h"
#include <vector>
#include <iostream>
using namespace std;

extern vector<Group> groupList;
extern int countGR, countDB;
extern int *arrDB, *arrGR, *cluster1, *cluster2;
extern FILE *joining1 , *joining2;
extern double EPS;
extern double EPS1;
extern int MINPOINTS;
extern int UNDEFINED;
extern int DIMENSION;    //value for the number of dimensions
extern int noOfPoints;
extern int GMINENTRIES;  
extern int GMAXENTRIES;
extern int GAUXMINENTRIES;
extern int GAUXMAXENTRIES; 
extern int RMINENTRIES;        //Minimum entries for node in RTree
extern int RMAXENTRIES;       //Maximum entries for node in RTree
extern double CELLSIZE;
extern int TIMESTAMP;
extern int CHOICE;
extern int CLUSTERID;
extern int IMMCOUNT;
extern int EPSCOUNT;
extern int BCELLID;
extern int GROUPID;
extern long long int DATAID;
extern double * MINGRIDSIZE;
extern double * MAXGRIDSIZE;
extern double *MINGRIDSIZEglobal, *MAXGRIDSIZEglobal;
extern double **MinGridBound, **MaxGridBound;
extern int numofThd;

extern int totalPoints;

CellDataHd initCellDataHd(){
    CellDataHd newCellDataHd = (CellDataHd)malloc(sizeof(struct celldatahd));
    ASSERT(newCellDataHd);
    newCellDataHd->count=0;
    newCellDataHd->first=NULL;
    newCellDataHd->isProcessed=FALSE;
    newCellDataHd->isDetermined=FALSE;
    return newCellDataHd;
}

Group initGroup(Data datapoint)
{
    Group newGroup = new group;
    //(Group) malloc(sizeof(struct group));

    newGroup->id= GROUPID++;
    newGroup->bcell = NULL;
    newGroup->master = datapoint->iData;
    newGroup->master_id = datapoint->id;
    newGroup->threshold = 0;
    newGroup->reachable_groups.clear();

    return newGroup;
}

BCell initBCell(Region rectTemp)
{
    // initializes an empty BCell
    BCell newBCell = (BCell) malloc(sizeof(*newBCell));
    newBCell->minOriginalBoundary = (Dimension)malloc(DIMENSION*sizeof(double));
    newBCell->maxOriginalBoundary = (Dimension)malloc(DIMENSION*sizeof(double));
    newBCell->minActualBoundary = (Dimension)malloc(DIMENSION*sizeof(double));
    newBCell->maxActualBoundary = (Dimension)malloc(DIMENSION*sizeof(double));
    newBCell->id= BCELLID++;

    int i =0;
    for(i = 0; i < DIMENSION; i++)
    {
        newBCell->minOriginalBoundary[i] = rectTemp->iBottomLeft[i];
        newBCell->maxOriginalBoundary[i] = rectTemp->iTopRight[i];
    }

    newBCell->cellDataList = initCellDataHd();

    return newBCell;

}

CellData initCellData(Data dataClstElem)
{
    CellData currCellData = (CellData) malloc(sizeof(*currCellData));
    ASSERT(currCellData);
    currCellData->data = dataClstElem;
    currCellData->data->core_tag = FALSE;

    currCellData->next = NULL;
    currCellData->nbhFlag=0;
   

    return currCellData;
}

GroupListHd initGroupListHd()
{
    groupList.clear();
    GroupListHd groupListHd = (GroupListHd) malloc(sizeof(*groupListHd));
    ASSERT(groupListHd);
    groupListHd->first = NULL;
    groupListHd->count=0;
    //printf("initialized succesfully %d\n", cellsListHd->count);
    return groupListHd;
}

BCellListHd initBCellListHd()
{
    BCellListHd cellsListHd = (BCellListHd) malloc(sizeof(*cellsListHd));
    ASSERT(cellsListHd);
    cellsListHd->first = NULL;
    cellsListHd->count=0;
    //printf("initialized succesfully %d\n", cellsListHd->count);
    return cellsListHd;
}

void freeCellsList(BCellListHd cellsListNbh)
{
    BCellListNode currListNode = cellsListNbh->first;
    BCellListNode nextListNode;
    
    while(currListNode!=NULL)
    {
        nextListNode = currListNode->next;
        free (currListNode);
        currListNode = nextListNode;
    }
    
    free(cellsListNbh);
    return;
}

void freeCellDataList(CellDataHd cellDataList)
{
    CellData currCellData = cellDataList->first;
    CellData nextCellData;
    while(currCellData!=NULL)
    {
        nextCellData=currCellData->next;
        free(currCellData);
        currCellData = nextCellData;
    }

    free(cellDataList);
    return;
}

void freeGroup(Group group)
{       
    group->reachable_groups.clear();
    free(group->master);

    free(group->bcell->minOriginalBoundary);
    free(group->bcell->maxOriginalBoundary);
    free(group->bcell->minActualBoundary);
    free(group->bcell->maxActualBoundary);

    freeAuxRTrees(group->bcell);
    CellData tempo = NULL;
    CellData temp = group->bcell->cellDataList->first;

    // while(temp != NULL)
    // {   
    //     tempo = temp->next;
    //     // if(temp->data->iData != NULL)free(temp->data->iData);
    //     // free(temp->data);
    //     // free(temp);
    //     temp = tempo;
    // }   

    free(group->bcell->cellDataList);
    free(group->bcell);

    free(group);
    return;
}

void freeGroupListHd(GroupListHd groupList)
{
    GroupListNode node = groupList->first;
    int i = 0;
    while(node != NULL)
    {
        GroupListNode temp = node;
        node = node->next;

        freeGroup(temp->groupElem);
        free(temp);
    }
    free(groupList);
    return;
}