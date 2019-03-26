#include "Grid.h"
#include <iostream>
using namespace std;

GHdrNd constructMainGRTree(DataHdr dataList1, GroupListHd cellsList, int gMinEntries, int gMaxEntries)
{
	GHdrNd GridRTree1 = (GHdrNd) populateGridRTree(dataList1,cellsList,gMinEntries,gMaxEntries);
	return GridRTree1;
}

void constructAuxRTrees(DataHdr dataList1, GroupListHd groupList)
{
	
	int cnt = 0,i;
	int j=0;
	Boolean cellFlag;

	GroupListNode currGroupNode = groupList->first;

	for(i=0;i<groupList->count;i++)
	{
		BCell currCell = currGroupNode->groupElem->bcell;			
		CellDataHd currCellDataList = currCell->cellDataList;
		CellData currCellData = currCellDataList->first;

		RHdrNd hdrNdTree = RinitHdrNd();
		hdrNdTree->ptrFirstNd = RinitLstNd(RinitIntNd(NULL, NULL));
		hdrNdTree->uiCnt++;
		hdrNdTree->ptrFirstNd->ptrChildLst = RinitHdrNd();
		
		while(currCellData!=NULL)
		{
			RinsertTree(hdrNdTree, RinitExtNd(currCellData->data));

				
			if(hdrNdTree->uiCnt > 1)
				hdrNdTree = RcreateRoot(hdrNdTree);

			currCellData=currCellData->next;			

		}

		currCell->auxCellRTree = hdrNdTree;

		currGroupNode = currGroupNode->next;
		
	}

	return;
}

void printAuxRTrees(BCellListHd cellsList)
{
	BCellListNode currCellNode=cellsList->first;
	BCell currBCell;
	int i =1;

	while(currCellNode!=NULL)
	{
		printf("\nPrinting %dth RTree\n",i++);

		currBCell = currCellNode->bCellElem;

		RprintTree(currBCell->auxCellRTree);

		currCellNode = currCellNode->next;
	}

	return;
}

void verifyAuxRTrees(BCellListHd cellsList)
{
	BCellListNode currCellNode=cellsList->first;
	BCell currBCell;
	int i =1;

	while(currCellNode!=NULL)
	{
		printf("\nVerifying %dth RTree\n",i++);

		currBCell = currCellNode->bCellElem;

		isCorrectRTree(currBCell->auxCellRTree);
		
		currCellNode = currCellNode->next;
	}

	return;
}

void freeAuxRTrees(BCell currCell)
{
	if(currCell->auxCellRTree !=NULL)
	{
		freeRTree(currCell->auxCellRTree);
	}	

	return;
}
