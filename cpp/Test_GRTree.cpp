#include "Test_GRTree.h"

int validateNeighborhood(RNbHdr n1,RNbHdr n2 )      
{
    //printf("R - %d\tGR - %d\n",n1->nbhCnt, n2->nbhCnt);
    if(n1->nbhCnt!=n2->nbhCnt)
    {
        printf("\tnbh count mismatch %d %d\n",n1->nbhCnt,n2->nbhCnt);
        getchar();
    }
    // int flag=0;
    // RNB nb_temp1 = n1->nbFirst,nb_next1,nb_temp2=n2->nbFirst,nb_next2;
    // while(nb_temp1!=NULL)
    // {
    //     nb_temp2=n2->nbFirst;
    //     flag=0;
    //     while(nb_temp2!=NULL)
    //     {
    //         if(nb_temp1->n == nb_temp2->n)
    //         {
    //             flag=1;
    //             break;
    //         }
    //         nb_temp2=nb_temp2->nbNext;
    //     }
    //     if(flag==0)
    //     {
    //         printf("point not found %d\n",nb_temp1->n );
    //         getchar();
    //     }
    //     //else
    //     //    printf("hi\n");

    //   nb_temp1=nb_temp1->nbNext;
    // }
    // nb_temp1=n2->nbFirst;
    // while(nb_temp1!=NULL)
    // {
    //     nb_temp2=n1->nbFirst;
    //     flag=0;
    //     while(nb_temp2!=NULL)
    //     {
    //         if(nb_temp1->n == nb_temp2->n)
    //         {
    //             flag=1;
    //             break;
    //         }
    //         nb_temp2=nb_temp2->nbNext;
    //     }
    //     if(flag==0)
    //     {
    //         printf("point not found %d\n",nb_temp1->n );
    //         getchar();
    //     }
    //     //else
    //     //    printf("hi\n");

    //   nb_temp1=nb_temp1->nbNext;
    // }

    return 0;
}

GHdrNd buildAuxGRTree(BCellListHd cellsList, int gMinEntries, int gMaxEntries)
{
    GHdrNd auxGRTree = populateAuxGridRTree(cellsList,gMinEntries,gMaxEntries);
    return auxGRTree;
}


void runTestCells(BCellListHd cellsList, GHdrNd GridRTree1, RHdrNd pointsRTree, DataHdr dataList1)
{
    // construct the GRTree
    // for every cell pick up every data point and initiate neighborhood queries
    // for all data points, do the neighborhood queries.

    BCellListNode currCellNode = cellsList->first;
    BCell currBCell;    
    BCellListHd cellsListNbh;
    Region cellRgnRect,epsExtRgnRect;
    GHdrNd auxGRTree;
    CellDataHd currCellDataList;
    CellData currCellData;
    Data currDataPoint;

    //Data dataTemp = (Data) malloc(sizeof(*dataTemp));

    int val; int count =1;
    int cellCount = 0;
    int cnt = 0;

    printf("Printing the No of Cells %d \n", cellsList->count);    
    //printNoOfCoreCells(cellsList);
    
    while (currCellNode!=NULL)
    {
        //printf("\nTesting for Cell %d\n",++cellCount);
        currBCell = currCellNode->bCellElem;
        // printf("printing current cell\n");
        // printCell(currBCell);

        cellRgnRect = createRegionofCell(currBCell);     
        epsExtRgnRect = getEpsExtendedRegion(cellRgnRect, (double)EPS);

        // printf("Printing Current Region \n");
        // GprintRegion(cellRgnRect);
        // printf("Printing Eps Extended Region \n");
        // GprintRegion(epsExtRgnRect);

        cellsListNbh = GgetCellsInRegion(GridRTree1,epsExtRgnRect,NULL);

        //printf("Printing Neighboring Cells:\n");
        //printCellsList(cellsListNbh);

        auxGRTree = buildAuxGRTree(cellsListNbh,GAUXMINENTRIES,GAUXMAXENTRIES);
        //printf("\nPrinting Aux GR Tree:\n");
        //GprintTree(auxGRTree);
        
        //RprintTree(pointsRTree);

        currCellDataList = currBCell->cellDataList;
        currCellData = currCellDataList->first;
        
        // printf("Printing neighbors of the current cell\n");
        //count = 1;
        while(currCellData!=NULL)
        {
            currDataPoint = currCellData->data;
            // printData(currDataPoint);
            // *dataTemp = *currDataPoint;
            GgetNeighborHood(auxGRTree, currDataPoint);

            // printf("\nPrinting Neighbors from aux\n\n");
            // printf("\n1 printing neighborhood of %lf %lf\n",currDataPoint->iData[0],currDataPoint->iData[1]);
            // RprintNbhLst(currDataPoint->neighbors, dataList1);
            // printf("\nPrinting Neighbors from original\n\n");
            // RgetNeighborHood(pointsRTree,dataTemp,0);
            // /printf("%d %d\n",currDataPoint->neighbors->nbhCnt, dataTemp->neighbors->nbhCnt );
            // printf("\n2 printing neighborhood of %lf %lf\n",dataTemp->iData[0],dataTemp->iData[1]);
            // RprintNbhLst(dataTemp->neighbors, dataList1);
           

            // validateNeighborhood(dataTemp->neighbors,currDataPoint->neighbors);
            RfreeNeighborhood(currDataPoint);
            // RfreeNeighborhood(dataTemp);
            // printf("printing neighbors of %dth datapoint\n", count++);
            // RprintNbhLst(currDataPoint->neighbors);
            // printf("\nend of printing neighbors\n");
            currCellData=currCellData->next;

            
        }
        // printf("Exiting Printing Neighborhood of the current cell");

        currCellNode=currCellNode->next;

        free(cellRgnRect->iBottomLeft);
        free(cellRgnRect->iTopRight);
        free(cellRgnRect);

        free(epsExtRgnRect->iBottomLeft);
        free(epsExtRgnRect->iTopRight);
        free(epsExtRgnRect);

        // free cellsListNbh

        freeCellsList(cellsListNbh);

        // free auxGRTree      

        freeGRTree(auxGRTree);

    }
        // cellRgnRect->iBottomLeft[0] = 50;
        // cellRgnRect->iBottomLeft[1] = 20;
        // cellRgnRect->iTopRight[0] = 60;
        // cellRgnRect->iTopRight[1] = 30;        
        // epsExtRgnRect = getEpsExtendedRegion(cellRgnRect, (double)EPS);        

        // printf("Printing Current Region \n");
        // GprintRegion(cellRgnRect);
        // printf("Printing Eps Extended Region \n");
        // GprintRegion(epsExtRgnRect);

        // cellsListNbh = GgetCellsInRegion(GridRTree1,epsExtRgnRect,NULL);

        // printf("Printing Neighboring Cells:\n");
        // printCellsList(cellsListNbh);        

        // auxGRTree = buildAuxGRTree(cellsListNbh,GAUXMINENTRIES,GAUXMAXENTRIES);
        // printf("\nPrinting Aux GR Tree:\n");
        // GprintTree(auxGRTree);
    
    // printf("\nPress Any Key to return from test function");
    // getchar();
    return;

}

void runTestPoints(BCellListHd cellsList, GHdrNd GridRTree1, RHdrNd pointsRTree, DataHdr dataList1)
{
    // construct the GRTree
    // for every cell pick up every data point and initiate neighborhood queries
    // for all data points, do the neighborhood queries.

    BCellListNode currCellNode = cellsList->first; 
    BCellListHd cellsListNbh;
    Region epsExtRgnRect;
    GHdrNd auxGRTree;
    CellDataHd currCellDataList;
    CellData currCellData;
    Data currDataPoint;

    

    // Data dataTemp = (Data) malloc(sizeof(*dataTemp));

    int val; int count =1;
    int cellCount = 0;
    int cnt = 0;

    Data temp = (Data) malloc(sizeof(*temp));
       

    //printf("Printing the No of Cells %d \n", cellsList->count);    
    //printNoOfCoreCells(cellsList);
    
    for(cnt=0;cnt<dataList1->uiCnt;cnt++)
    {
        //printf("Printing for %dth point\n",cnt+1);
        currDataPoint = dataList1->dataClstElem+cnt;
        *temp = *currDataPoint; 
        epsExtRgnRect = getEpsExtendedRegionPoint(currDataPoint,EPS);
        cellsListNbh = GgetCellsInRegion(GridRTree1,epsExtRgnRect,NULL);
        
        currCellNode = cellsListNbh->first;   
        currDataPoint->neighbors = RinitNbHdr();    

        while(currCellNode!=NULL)
        {
            val = RgetNeighborHood(currCellNode->bCellElem->auxCellRTree,temp,0);

            //nbh of currDataPoint is to be appended to the temp
            appendNbh(currDataPoint,temp);
            currCellNode = currCellNode->next;
        }

        //free(temp);

        RfreeNeighborhood(currDataPoint);

        free(epsExtRgnRect->iBottomLeft);
        free(epsExtRgnRect->iTopRight);
        free(epsExtRgnRect);
        freeCellsList(cellsListNbh);
    }

    free(temp);
   
    return;

}
