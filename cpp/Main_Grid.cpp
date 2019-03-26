#include "Main_Grid.h"

GHdrNd constructMainGRTree(DataHdr dataList1, BCellListHd cellsList, int gMinEntries, int gMaxEntries)
{
	GHdrNd GridRTree1 = (GHdrNd) populateGridRTree(dataList1,cellsList,gMinEntries,gMaxEntries);
	return GridRTree1;
}

void constructAuxRTrees(DataHdr dataList1, BCellListHd cellsList)
{
	
	int cnt = 0,i;
	int j=0;
	Boolean cellFlag;
    //retrieve element one by one and insert them into tree invoking create root incase of root split
	
	BCellListNode currCellNode = cellsList->first;

	for(i=0;i<cellsList->count;i++)
	{
		BCell currCell = currCellNode->bCellElem;			
		CellDataHd currCellDataList = currCell->cellDataList;
		CellData currCellData = currCellDataList->first;

		RHdrNd hdrNdTree = RinitHdrNd();
		hdrNdTree->ptrFirstNd = RinitLstNd(RinitIntNd(NULL, NULL));
		hdrNdTree->uiCnt++;
		hdrNdTree->ptrFirstNd->ptrChildLst = RinitHdrNd();
		
		while(currCellData!=NULL)
		{
			/*printf("inserting %d %lf %lf\n",currCellData->data->iNum, currCellData->data->iData[0], currCellData->data->iData[1] );
			if(currCellData->data->iNum==58 || currCellData->data->iNum==62 || currCellData->data->iNum==63)
				getchar(); */
			
			RinsertTree(hdrNdTree, RinitExtNd(currCellData->data));
			//printf("root count = %d\n",hdrNdTree->uiCnt );
			// if(currCellData->data->iNum==58 || currCellData->data->iNum==62 || currCellData->data->iNum==63)
			// {
			// 	getchar();
			// 	RprintTree(hdrNdTree);
			// }
				
			if(hdrNdTree->uiCnt > 1)
				hdrNdTree = RcreateRoot(hdrNdTree);

			// if(currCellData->data->iNum==58 || currCellData->data->iNum==62 || currCellData->data->iNum==63)
			// {
			// 	getchar();
			// 	RprintTree(hdrNdTree);
			// }

			currCellData=currCellData->next;			

		}

		currCell->auxCellRTree = hdrNdTree;

		currCellNode = currCellNode->next;
		
	}

	return;
}

void constructAuxRTreesWrtDense(DataHdr dataList1, BCellListHd cellsList, int denseCount)
{
	
	int cnt = 0,i;
	int j=0;
	Boolean cellFlag;
    //retrieve element one by one and insert them into tree invoking create root incase of root split
	
	BCellListNode currCellNode = cellsList->first;

	for(i=0;i<cellsList->count;i++)
	{
		BCell currCell = currCellNode->bCellElem;
		if(currCell->cellDataList->count > denseCount)
		{
			currCell->isDenseCell = TRUE;

			BCellListHd currCellsList = initBCellListHd();

			currCell->auxCellGRTree = populateCellGridRTree(currCell->cellDataList, currCellsList, currCell, GAUXMINENTRIES, GAUXMAXENTRIES);
			currCell->auxCellsList = currCellsList;

			constructAuxRTrees(dataList1,currCellsList);
		}	
		else
		{
			currCell->isDenseCell = FALSE;

			CellDataHd currCellDataList = currCell->cellDataList;
			CellData currCellData = currCellDataList->first;

			RHdrNd hdrNdTree = RinitHdrNd();
			hdrNdTree->ptrFirstNd = RinitLstNd(RinitIntNd(NULL, NULL));
			hdrNdTree->uiCnt++;
			hdrNdTree->ptrFirstNd->ptrChildLst = RinitHdrNd();
		
			while(currCellData!=NULL)
			{
				/*printf("inserting %d %lf %lf\n",currCellData->data->iNum, currCellData->data->iData[0], currCellData->data->iData[1] );
				if(currCellData->data->iNum==58 || currCellData->data->iNum==62 || currCellData->data->iNum==63)
					getchar(); */
				
				RinsertTree(hdrNdTree, RinitExtNd(currCellData->data));
				//printf("root count = %d\n",hdrNdTree->uiCnt );
				// if(currCellData->data->iNum==58 || currCellData->data->iNum==62 || currCellData->data->iNum==63)
				// {
				// 	getchar();
				// 	RprintTree(hdrNdTree);
				// }
					
				if(hdrNdTree->uiCnt > 1)
					hdrNdTree = RcreateRoot(hdrNdTree);

				// if(currCellData->data->iNum==58 || currCellData->data->iNum==62 || currCellData->data->iNum==63)
				// {
				// 	getchar();
				// 	RprintTree(hdrNdTree);
				// }

				currCellData=currCellData->next;			

			}

			currCell->auxCellRTree = hdrNdTree;

		}				

		currCellNode = currCellNode->next;
		
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

void freeAuxRTrees(BCellListHd cellsList)
{
	
	int cnt = 0,i;
	int j=0;
	Boolean cellFlag; int rmin, rmax;
    //retrieve element one by one and insert them into tree invoking create root incase of root split
	
	BCellListNode currCellNode = cellsList->first;
	BCell currCell;

	for(i=0;i<cellsList->count;i++)
	{
		currCell = currCellNode->bCellElem;

		if(currCell->auxCellRTree !=NULL)
		{
			freeRTree(currCell->auxCellRTree);			
		}

		if(currCell->auxCellGRTree !=NULL)	
		{
			freeGRTree(currCell->auxCellGRTree);
			freeAllBCells(currCell->auxCellsList);
		}

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

void CellsPointsTest(BCellListHd cellsList, int limit)
{
	BCellListNode currCellNode=cellsList->first;
	BCell currBCell;
	int i =1;

	int sparse = 0; int dense =0;

	while(currCellNode!=NULL)
	{
		currBCell = currCellNode->bCellElem;

		if (currBCell->cellDataList->count > limit)
		{
			dense++;
		}
		else
		{
			sparse++;
		}

	
		currCellNode = currCellNode->next;
	}

	printf("No of Dense Cells: %d\n", dense);
	printf("No of Sparse Cells: %d\n", sparse);

	return;
}

void captureDenseCells(BCellListHd cellsList, BCellListHd denseCellsList, int limit)
{
	BCellListNode currCellNode=cellsList->first;
	BCell currBCell;
	int i =1;

	int sparse = 0; int dense =0;

	while(currCellNode!=NULL)
	{
		currBCell = currCellNode->bCellElem;

		if (currBCell->cellDataList->count >= limit)
		{
			insertBCellIntoBCellList(currBCell,denseCellsList);
		}
	
		currCellNode = currCellNode->next;
	}

	return;

}

DataHdr createDataListFromCellList(BCellListHd cellsList)
{
	int noOfPoints=0; int i=0;
	BCell currBCell;

	BCellListNode currCellNode=cellsList->first;
	while(currCellNode!=NULL)
	{
		currBCell = currCellNode->bCellElem;
		noOfPoints = noOfPoints + currBCell->cellDataList->count;			
		currCellNode = currCellNode->next;
	}

	DataHdr tempDataList = initDataHdr(noOfPoints);
	DataPoint dataPointTemp = NULL;
	CellData currCellData;
	int iDim=0;

	currCellNode=cellsList->first;
	while(currCellNode!=NULL)
	{
		
		currBCell = currCellNode->bCellElem;
		currCellData =  currBCell->cellDataList->first;
		
		while(currCellData!=NULL)
		{
			dataPointTemp = (DataPoint)malloc(DIMENSION *sizeof(dataPoint));

        	for(iDim = 0; iDim < DIMENSION; iDim++)
        		dataPointTemp[iDim] = currCellData->data->iData[iDim];
        
        	insertDataLstElem(tempDataList, dataPointTemp);
			currCellData = currCellData->next;
		}

		currCellNode = currCellNode->next;
	}

	return tempDataList;

}

BCellListHd captureDenseCellsIncludingNeighbors(GHdrNd GridRTree1,BCellListHd cellsList,BCellListHd denseCellsList)
{
	BCellListNode currBCellListNode = denseCellsList->first;
	BCellListNode tempBCellListNode;
	BCell currBCell;
	CellData currCellData;
	int counter=0;
	int i;
	double randCoords[DIMENSION];

	BCellListHd nbhCellsList = initBCellListHd();
	
	while(currBCellListNode!=NULL)
	{
		Region cellRgn = createRegionofCell(currBCellListNode->bCellElem);
		cellRgn=getEpsExtendedRegion(cellRgn,CELLSIZE);
		BCellListHd cellNbhCells = GgetCellsInRegion(GridRTree1,cellRgn,NULL);

		tempBCellListNode = cellNbhCells->first;;

		while(tempBCellListNode!=NULL)
		{

			currBCell = tempBCellListNode->bCellElem;

			if(currBCell->cellDataList->count < 10000)
			{
				counter = abs(10000 - currBCell->cellDataList->count);

	 			while(counter > 0)
	 			{
	 				for(i = 0; i<DIMENSION; i++)
		 			{
		 				randCoords[i] = fRand(currBCell->minOriginalBoundary[i], currBCell->maxOriginalBoundary[i]);
		 			}

		 			Data currData = (Data) malloc(sizeof(*currData));
		 			currData->iData = (DataPoint) malloc(DIMENSION *sizeof(dataPoint));
		 			for(i=0;i<DIMENSION;i++)
		 			{
		 				currData->iData[i] = randCoords[i];
		 			}

		 			insertPointintoBCell(currBCell,currData);
		 			counter--;
	 			}



				// currCellData = currBCell->cellDataList->first;

				// while(currBCell->cellDataList->count > 20)
	 		// 	{
	 		// 		currCellData = currCellData->next;
	 		// 		currBCell->cellDataList->count--;
	 		// 	}

	 			//currBCell->cellDataList->first = currCellData;

			}
			//insertBCellIntoBCellList(tempBCellListNode->bCellElem,nbhCellsList);
			tempBCellListNode = tempBCellListNode->next;
		}

		currBCellListNode = currBCellListNode->next;
	}	

	return cellsList;

}


BCellListHd removeDuplicateCells(BCellListHd denseCellsList)
{
	BCellListNode currBCellListNode = denseCellsList->first;
	BCell currBCell;
	BCellListNode tempBCellListNode;
	BCellListNode pred;
	BCellListNode succ;
	Boolean flag;

	BCellListHd newList = initBCellListHd();

	currBCell = currBCellListNode->bCellElem;
	insertBCellIntoBCellList(currBCell,newList);

	currBCellListNode = currBCellListNode->next;

	while(currBCellListNode!=NULL)
	{
		currBCell = currBCellListNode->bCellElem;

		tempBCellListNode = newList->first;

		while(tempBCellListNode != NULL)
		{
			if((isSameCells(currBCell,tempBCellListNode->bCellElem)))
			{
				flag = TRUE;
				break;
			}
			tempBCellListNode = tempBCellListNode->next;
		}

		if(flag == FALSE)
		{
			insertBCellIntoBCellList(currBCell,newList);
		}

		currBCellListNode = currBCellListNode->next;
	}

	return newList;
}

Boolean isSameCells(BCell first, BCell second)
{
	int i = 0 ; int counter = 0;

	for(i = 0; i < DIMENSION;i++)
	{
		if(first->minActualBoundary[i] == second->minActualBoundary[i])
		{
			counter++;
		}
		if(first->maxActualBoundary[i] == second->maxActualBoundary[i])
		{
			counter++;
		}

	}

	if(counter == (DIMENSION  * 2))
	{
		return TRUE;
	}

	return FALSE;

}


void increasePointsToInDenseCells(BCellListHd denseCellsList,int limit)
{
	BCellListNode currBCellListNode = denseCellsList->first;
	BCell currBCell;
 	BCellListNode tempBCellListNode;
 	Data currData;
 	DataPoint currDataPoint;

 	int i =0;

 	srand((unsigned)time(NULL));

 	int counter=0;

 	double randCoords[DIMENSION];
 	
 	while(currBCellListNode!=NULL)
 	{	
		currBCell = currBCellListNode->bCellElem;

 		if(currBCell->cellDataList->count<=limit)
 		{
 			counter = limit - currBCell->cellDataList->count;

 			while(counter > 0)
 			{
 				for(i = 0; i<DIMENSION; i++)
	 			{
	 				randCoords[i] = fRand(currBCell->minOriginalBoundary[i], currBCell->maxOriginalBoundary[i]);
	 			}

	 			Data currData = (Data) malloc(sizeof(*currData));
	 			currData->iData = (DataPoint) malloc(DIMENSION *sizeof(dataPoint));
	 			for(i=0;i<DIMENSION;i++)
	 			{
	 				currData->iData[i] = randCoords[i];
	 			}

	 			insertPointintoBCell(currBCell,currData);
	 			counter--;
 			}

 		}

 		currBCellListNode = currBCellListNode->next;
 		counter = 0;

 	}

}

void reducePointsOfSparseNeighbors(BCellListHd denseCellsList, int limit)
{
	BCellListNode currBCellListNode = denseCellsList->first;
	BCell currBCell;
 	BCellListNode tempBCellListNode;
 	Data currData;
 	DataPoint currDataPoint;

 	int i =0;

 	srand((unsigned)time(NULL));

 	int counter=0;

 	double randCoords[DIMENSION];

 	//CellDataHd currCellDatalist;
 	CellData currCellData;

 	
 	while(currBCellListNode!=NULL)
 	{	
		currBCell = currBCellListNode->bCellElem;

 		if(currBCell->cellDataList->count>limit)
 		{
 			counter = abs(currBCell->cellDataList->count - limit);

 			//currCellDatalist = currBCell->cellDataList;
 			currCellData = currBCell->cellDataList->first;

 			// remove points from the cellDataList

 			while(currBCell->cellDataList->count > 20)
 			{
 				currCellData = currCellData->next;
 				currBCell->cellDataList->count--;
 			}

 			//currBCell->cellDataList->count = 20;
 			currBCell->cellDataList->first = currCellData;

 		}

 		currBCellListNode = currBCellListNode->next;
 		counter = 0;

 	}

}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void writeSelectedPointsToFileFromCellsList(BCellListHd denseCellsList)
{
	BCellListNode currBCellListNode = denseCellsList->first;

	BCell currBCell; int j=0;

	FILE *filP = fopen("selectedPoints", "w");

	fprintf(filP, "%d\n", denseCellsList->count);
    fprintf(filP, "%d\n", DIMENSION);


	while(currBCellListNode!=NULL)
	{
		currBCell = currBCellListNode->bCellElem;

		for (j = 0; j < DIMENSION; j++)
	    {
            fprintf(filP, "%lf", currBCell->cellDataList->first->data->iData[j]);
            
            if(j!=(DIMENSION-1))
            {
                fprintf(filP," ");
            }
        }

	    fprintf(filP,"\n");
	    
		currBCellListNode = currBCellListNode->next;
	}

	fclose(filP);

	return;
}

BCellListHd captureSelectedCells(BCellListHd cellsList, BCellListHd denseCellsList, DataHdr sampledDataList, GHdrNd GridRTree1)
{
	int i = 0;

	Region currRegion;
	BCell currBCell;

	for(i=0;i<sampledDataList->uiCnt;i++)
	{
		currRegion = createCellRegOfPoint(sampledDataList->dataClstElem + i);
		currBCell = GfindCell(GridRTree1,currRegion);

		insertBCellIntoBCellList(currBCell,denseCellsList);

	}

	return denseCellsList;

}

BCellListHd segregateSparseCellsWhoseNBHIsSparse(BCellListHd cellsList, GHdrNd GridRTree1)
{
	BCellListHd sparseCellsList = initBCellListHd();

	BCellListNode currBCellListNode = cellsList->first;

	BCellListHd tempCellsList;

	BCellListNode tempBCellListNode;

	Region tempRegion; Region tempRegion2;

	int globalCount = 0; Boolean flag;

	while(currBCellListNode!=NULL)
	{
		tempRegion = createRegionofCell(currBCellListNode->bCellElem);
		tempRegion2 = getEpsExtendedRegion(tempRegion,CELLSIZE);
		
		tempCellsList = GgetCellsInRegion(GridRTree1, tempRegion2, NULL);

		tempBCellListNode = tempCellsList->first;

		while(tempBCellListNode!=NULL)
		{
			if(tempBCellListNode->bCellElem->cellDataList->count >= 3000)
			{
				flag=TRUE;
				break;
			}

			tempBCellListNode = tempBCellListNode->next;
		}

		if(flag == FALSE)
		{
			insertBCellIntoBCellList(currBCellListNode->bCellElem,sparseCellsList);
			globalCount = globalCount + currBCellListNode->bCellElem->cellDataList->count;
		}

		if(globalCount >= 80000)
		{
			break;
		}

		currBCellListNode = currBCellListNode->next;
		flag = FALSE;

	}

	return sparseCellsList;

}

int main(int argc, char **argv)
{
	//DIMENSION = 2;

// ./output $path1$1 output_$1\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS.txt $GMINENTRIES $GMAXENTRIES $RMINENTRIES $RMAXENTRIES $CELLSIZE $Epsilon $MINPOINTS $UNDEFINED

	EPS = atof(argv[10]);
	//printf("EPS is %lf\n", EPS);
	GMINENTRIES = atoi(argv[3]);
	//printf("GMINENTRIES is %d\n", GMINENTRIES);
	GMAXENTRIES = atoi(argv[4]);
	//printf("GMAXENTRIES is %d\n", GMAXENTRIES);
	GAUXMINENTRIES = atoi(argv[5]);
	//printf("GAUXMINENTRIES is %d\n", GAUXMINENTRIES);
	GAUXMAXENTRIES = atoi(argv[6]);
	//printf("GAUXMAXENTRIES is %d\n", GAUXMAXENTRIES);
	RMINENTRIES = atoi(argv[7]);
	//printf("RMINENTRIES is %d\n", RMINENTRIES);
	RMAXENTRIES = atoi(argv[8]);
	//printf("RMAXENTRIES is %d\n", RMAXENTRIES);
	MINPOINTS = atoi(argv[11]);
	//printf("MINPOINTS is %d\n", MINPOINTS);
	CELLSIZE = atof(argv[9]);
	//printf("CELLSIZE is %lf\n", CELLSIZE);
	UNDEFINED = atoi(argv[12]);
	//printf("UNDEFINED is %d\n", UNDEFINED); 
	int k = atoi(argv[13]);
	// printf("k is %d\n", k); 

	char * QT = (char *)argv[15];
	//printf("%s\n",argv[15]);

/*	EPS = 7;
	GMINENTRIES = 4;
	GMAXENTRIES = 8;
	RMINENTRIES = 2;
	RMAXENTRIES = 4;
	MINPOINTS = 3;
	CELLSIZE = 10;
	UNDEFINED = 100000000; */
	


	// declarations for dataList and RTree
	DataHdr dataList1;
	DataHdr tempDataList;
	GHdrNd GridRTree1;
	DataHdr sampledDataList;

	int noOfPoints,i,iCnt;
	dataList1 = readData(argv[1]); //passing filename
	//printf("Reading1 Successfull\n");
	//printDataLst(dataList1);

	sampledDataList = readSampledData(argv[14]);
	//printf("Reading2 Successfull\n");

	//printMinGridSize();
	//printMaxGridSize();
	//printf("hello 1\n");
	// initialize a new cells list
	BCellListHd cellsList = initBCellListHd();
	//BCellListHd sparseCellsList = initBCellListHd();
	BCellListHd denseCellsList = initBCellListHd();
	
	RHdrNd testRTree;
	//testRTree = RbuildRTree(dataList1);

	//RprintTree(testRTree);
	//printf("hello 2\n");
	//populate grid based R Tree and the cellsList
	GridRTree1 = constructMainGRTree(dataList1,cellsList,GMINENTRIES,GMAXENTRIES);
	//isCorrectGRTree(GridRTree1);
	//printf("\nVerified Main GridRTree\n");

	//CellsPointsTest(cellsList,3000);

	
	//tempDataList = createDataListFromCellList(denseCellsList);
	//writeDataListToFile(tempDataList, "3D_spatial_network_transformed_dense_cells_more_than_3000pts");
	//GprintTree(GridRTree1);
	//getchar();
	
	//printf("hello 3\n");
	//constructAuxRTrees(dataList1,cellsList);
	constructAuxRTreesWrtDense(dataList1,cellsList,3500);
	//verifyAuxRTrees(cellsList);
	//printf("\nVerified Aux RTrees\n");
	//printAuxRTrees(cellsList);
	//printf("hello 4\n");
	//getchar();
	
	
	denseCellsList=captureSelectedCells(cellsList,denseCellsList,sampledDataList, GridRTree1);
	// captureDenseCells(cellsList, denseCellsList, 3000);
	// cellsList = captureDenseCellsIncludingNeighbors(GridRTree1,cellsList,denseCellsList);	
	// increasePointsToInDenseCells(denseCellsList,10000);
	// reducePointsOfSparseNeighbors(denseCellsList,20);
	// writeSelectedPointsToFileFromCellsList(denseCellsList);
	// BCellListHd sparseCellsList = segregateSparseCellsWhoseNBHIsSparse(cellsList,GridRTree1);
	// tempDataList = createDataListFromCellList(sparseCellsList);
	// writeDataListToFile(tempDataList, "deluciaD32lac_sparse_cells_less_than_3000pts_after_increasing_densecells_and_nbhs_to_10000");
	// writeDataListToFile(tempDataList, "3D_spatial_network_transformed_dense_cells_greater_than_3000pts_increased_to_1000pts_withoutNBH");
 
	// runTestOneNN(dataList1,GridRTree1);
	// runTestKNN(dataList1,GridRTree1,4);
	// runTestKNNOnSample(dataList1,sampledDataList,GridRTree1,k);

	// runTestCells(cellsList,GridRTree1,testRTree,sampledDataList);

	if(strcmp(QT,"CELLWISENBH") == 0)
	{
		//printf("\nGoing for CELLWISENBH\n");
		runTestCellsWrtDense(denseCellsList,GridRTree1,testRTree,dataList1);
		//freeCellsList(denseCellsList);
	}
	else if(strcmp(QT,"POINTWISENBH") == 0)
	{
		//printf("\nGoing for POINTWISENBH\n");
		//runTestPoints(cellsList,GridRTree1,testRTree,sampledDataList);
		runTestPointsWrtDense(cellsList,GridRTree1,testRTree,sampledDataList);

	}
	else if(strcmp(QT,"KNN") == 0)
	{
		//printf("\nGoing for KNN\n");
		//runTestKNNOnSample(dataList1,sampledDataList,GridRTree1,k);
		runTestKNNOnSampleWrtDense(dataList1,sampledDataList,GridRTree1,k);
	}

	// printf("\nPress Any Key to exit");
	// getchar();

	//freeAuxRTrees(cellsList);
	//freeGRTree(GridRTree1);
	//freeAllBCells(cellsList);
	//freeCellsList(denseCellsList);
	//freeDataList(dataList1);
	//freeDataList(sampledDataList);

	//free(MINGRIDSIZE);
	//free(MAXGRIDSIZE);

    return 0;
}
