#include "RTree.h"
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

extern vector<Group> groupList;

extern int totalPoints;

// extern Cluster cluster;

Region RinitRgnRect(Dimension iBottomLeft, Dimension iTopRight)
{   //initializes the rectangle with the given bottom left and top right corners
    //if the values for the corners are specified NULL, initializes a rectangle with origin as co-ordinates for both corners.
	Region rgnRect = (Region)malloc(sizeof(*rgnRect));
    if(rgnRect == NULL)
		return NULL;
		
    if(iBottomLeft != NULL)
		rgnRect->iBottomLeft = iBottomLeft;
	else
    //rgnRect->iBottomLeft = (Dimension) calloc( DIMENSION, sizeof( dimension ) );
		rgnRect->iBottomLeft = (Dimension)malloc(sizeof(double)*DIMENSION);

	if(rgnRect->iBottomLeft == NULL)
    {   free(rgnRect);
		return NULL;
    }

	if(iTopRight != NULL)
		rgnRect->iTopRight = iTopRight;
	else
		rgnRect->iTopRight = (Dimension)malloc(sizeof(double)*DIMENSION);

	if(rgnRect->iTopRight == NULL)
    {   
    	if(rgnRect->iBottomLeft == NULL);
		free(rgnRect);
		return NULL;
	}

	return rgnRect;
}

void RsetRect(RLstNd lstNd, RTreeNode tnInfo)
{   // copies the data in the tree node tnInfo to lstNd
    int iCnt = 0;
    switch(tnInfo->ndType)
    {   case INTNODE:
		//incase of internal node copy the bottom left and top right corners
		for(iCnt = 0; iCnt < DIMENSION; iCnt++)
        {   lstNd->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] = tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt];
			lstNd->tnInfo->tdInfo->rgnRect->iTopRight[iCnt] = tnInfo->tdInfo->rgnRect->iTopRight[iCnt];
     	}
		break;
        case EXTNODE:
		// in case of external node copy the data element
		for(iCnt = 0; iCnt < DIMENSION; iCnt++)
        {   lstNd->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] = tnInfo->tdInfo->dataClstElem->iData[iCnt];
			lstNd->tnInfo->tdInfo->rgnRect->iTopRight[iCnt] = tnInfo->tdInfo->dataClstElem->iData[iCnt];
		}
		break;
	}

	return;
}

RTreeNode RinitIntNd(Dimension iBottomLeft, Dimension iTopRight)
{   //intializes internal node of a Tree with rectangle whose bottom left and topright corners are given

	Region rgnRect = RinitRgnRect(iBottomLeft, iTopRight);
	// initializes a rectangle with the given coordonates for the bottom left and top right corners

	RTreeData tdInfo = (RTreeData)malloc(sizeof(*tdInfo));

	if(tdInfo == NULL)
		return NULL;

	tdInfo->rgnRect = rgnRect;
    RTreeNode tnInfo = (RTreeNode)malloc(sizeof(*tnInfo));

	if(tnInfo == NULL)
		return NULL;

	tnInfo->ndType = INTNODE;
	tnInfo->tdInfo = tdInfo;

	return tnInfo;
}

RTreeNode RinitExtNd(Data dataClstElem)
{  
   // initializes the external node of a tree with Data
   if(dataClstElem == NULL)
		return NULL;

	RTreeNode tnInfo = (RTreeNode)malloc(sizeof(*tnInfo));

	if(tnInfo == NULL)
		return NULL;

	RTreeData tdInfo = (RTreeData)malloc(sizeof(*tdInfo));
    if(tdInfo == NULL)
    {   free(tnInfo);
		return NULL;
	}

	tdInfo->dataClstElem = dataClstElem;	// Data
    tnInfo->ndType = EXTNODE;	// external node
	tnInfo->tdInfo = tdInfo;

	return tnInfo;
}

RHdrNd RcreateRoot(RHdrNd hdrNdTree)
{   //in case of root split this is called to create a new root
    RHdrNd hdrNdRoot = RinitHdrNd();
    Dimension iBottomLeft = (Dimension)calloc(DIMENSION, sizeof(dimension));
	Dimension iTopRight = (Dimension)calloc(DIMENSION,sizeof(dimension));

	RLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
	int iCnt = 0;
	Boolean bIsFirst = TRUE;

   //set the bottom left and top right corners for the new root
	while(lstNdTemp != NULL)
    {	for(iCnt = 0; iCnt < DIMENSION; iCnt++)
        {   if(bIsFirst)
            {   iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt];
				iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt];
         	}
			else
            {   if(lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] < iBottomLeft[iCnt])
					iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt];
				if(lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt] > iTopRight[iCnt])
					iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt];
    	    }
		}
		lstNdTemp = lstNdTemp->ptrNextNd;
		bIsFirst = FALSE;
	}

//initialize a node with the bottomleft and top right corners obtained and insert into the list
	hdrNdRoot->ptrFirstNd = RinitLstNd(RinitIntNd(iBottomLeft, iTopRight));
    hdrNdRoot->ptrFirstNd->ptrChildLst = hdrNdTree;
	hdrNdRoot->uiCnt = 1;

	return hdrNdRoot;
}

Boolean RexpansionArea(Region rgnRect, RTreeNode tnInfo, Double ptrDMinExp, Region rgnNewRect)
{   // calculates if the area by which the rgnRect should be enlarged so as to include the tnInfo is less than the value pointed by ptrDMinExp and return TRUE and assigns rgnNewRect with the new enlarged rectangle.
    int iCnt = 0;
    Region rgnRectTemp = RinitRgnRect(NULL, NULL);
    for(iCnt = 0; iCnt < DIMENSION; iCnt++)
    {   switch(tnInfo->ndType)
        {   case INTNODE:
            //assign least of bottom left corner along each dimension to rgnRectTemp
			rgnRectTemp->iTopRight[iCnt] = (tnInfo->tdInfo->rgnRect->iTopRight[iCnt] > rgnRect->iTopRight[iCnt]) ? tnInfo->tdInfo->rgnRect->iTopRight[iCnt] : rgnRect->iTopRight[iCnt];
			rgnRectTemp->iBottomLeft[iCnt] = (tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] < rgnRect->iBottomLeft[iCnt]) ? tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] : rgnRect->iBottomLeft[iCnt];
            break;

		    case EXTNODE:
           //assign maximum of top right corner along each dimension to rgnRectTemp
			rgnRectTemp->iTopRight[iCnt] = (tnInfo->tdInfo->dataClstElem->iData[iCnt] > rgnRect->iTopRight[iCnt]) ? tnInfo->tdInfo->dataClstElem->iData[iCnt] : rgnRect->iTopRight[iCnt];
			rgnRectTemp->iBottomLeft[iCnt] = (tnInfo->tdInfo->dataClstElem->iData[iCnt] < rgnRect->iBottomLeft[iCnt]) ? tnInfo->tdInfo->dataClstElem->iData[iCnt] : rgnRect->iBottomLeft[iCnt];
			break;
		}
	}
    //calculate the difference in area for new rectangle and old rectangle
	double dExp = Rarea(rgnRectTemp) - Rarea(rgnRect);
    //in case there no value to compare ( -1 ) or incase the value is less than the value to be comparedcopy the rgnRectTemp to rgnRectNew to Return it.
	if(*ptrDMinExp == -1 || dExp <= *ptrDMinExp)
    {   
    	if(dExp == *ptrDMinExp)
			*ptrDMinExp = 0 - dExp;
		else
			*ptrDMinExp = dExp;
        for(iCnt =0; iCnt< DIMENSION; iCnt++)
        {	rgnNewRect->iBottomLeft[iCnt] = rgnRectTemp->iBottomLeft[iCnt];
			rgnNewRect->iTopRight[iCnt] = rgnRectTemp->iTopRight[iCnt];
		}
		free(rgnRectTemp->iBottomLeft);
		free(rgnRectTemp->iTopRight);
		free(rgnRectTemp);
    //area to be enlarged is less than the previous value
		return TRUE;
	}

	free(rgnRectTemp->iBottomLeft);
	free(rgnRectTemp->iTopRight);
	free(rgnRectTemp);
    //area to be enlarged is not less than the previous value
	return FALSE;
}

double Rarea(Region rgnRect)
{   //calcluates the area of rectangle
    if(rgnRect == NULL)
		return 0;
    double dArea = 1;
	int iCnt = 0;
    //multiply values along each dimension
	for(iCnt = 0; iCnt < DIMENSION; iCnt++)
		dArea = dArea * (rgnRect->iTopRight[iCnt] - rgnRect->iBottomLeft[iCnt]);
	return dArea;
}

RLstNd RpickChild(RHdrNd ptrChildLst, RTreeNode tnInfo)
{  // decides which node among the child nodes to be picked for insertion and returns a pointer to that node
    if(ptrChildLst == NULL)
		return NULL;

	RLstNd lstNdTemp = ptrChildLst->ptrFirstNd;
	//GprintRegion( lstNdTemp->tnInfo->tdInfo->rgnRect );
	RLstNd lstNdChild = NULL;
	double dMinExp = -1;
	int iCnt;
    Region rgnNewRect = RinitRgnRect(NULL, NULL);
    Region rgnFinalRect = RinitRgnRect(NULL, NULL);
    // for each child child in the list of child nodes do the following
	while(lstNdTemp != NULL)
    {   //call the expansionArea function to determine the are by which the child node has to enlarged to accomodate the new point or region.
		if(RexpansionArea(lstNdTemp->tnInfo->tdInfo->rgnRect, tnInfo, &dMinExp, rgnNewRect))
        {//if the expansionArea return true mark the node to be the one that might be picked. if the expansion is same as one of the previous nodes then compare the ares of the current noe and the previous node and pick the one with least area.
               if(dMinExp < 0)
               {     dMinExp = 0 - dMinExp;
                     if(Rarea(lstNdChild->tnInfo->tdInfo->rgnRect) > Rarea(lstNdTemp->tnInfo->tdInfo->rgnRect))
                     {
                     	     lstNdChild = lstNdTemp;
                     	     for(iCnt =0; iCnt< DIMENSION; iCnt++)
					         {	rgnFinalRect->iBottomLeft[iCnt] = rgnNewRect->iBottomLeft[iCnt];
								rgnFinalRect->iTopRight[iCnt] = rgnNewRect->iTopRight[iCnt];
							 }
                     }
               }
			   else
			   {
				   lstNdChild = lstNdTemp;
				   for(iCnt =0; iCnt< DIMENSION; iCnt++)
		         	{	rgnFinalRect->iBottomLeft[iCnt] = rgnNewRect->iBottomLeft[iCnt];
						rgnFinalRect->iTopRight[iCnt] = rgnNewRect->iTopRight[iCnt];
				 	}
				}
				//lstNdChild = lstNdTemp;
		}
    	lstNdTemp = lstNdTemp->ptrNextNd;
	}
    //for the node that is picked assign the region pointed by new rectangle region and return the node
   //GprintRegion( lstNdChild->tnInfo->tdInfo->rgnRect );
    Region rgnRectTemp = lstNdChild->tnInfo->tdInfo->rgnRect;
	lstNdChild->tnInfo->tdInfo->rgnRect = rgnFinalRect;
	//GprintRegion( lstNdChild->tnInfo->tdInfo->rgnRect );
	//rgnNewRect = NULL;

	
	free(rgnRectTemp->iBottomLeft);
	free(rgnRectTemp->iTopRight);
	free(rgnRectTemp);

	free(rgnNewRect->iBottomLeft);
	free(rgnNewRect->iTopRight);
	free(rgnNewRect);

	return lstNdChild;
}

void RpickSeeds(RHdrNd ptrChildLst, RLstNd *lstNdChildOne, RLstNd *lstNdChildTwo)
{    //in case a node has to be split pick the two child nodes that are used to create new child lists
	if(ptrChildLst == NULL)
		return;
    //used to store pointers to nodes of maximum bottom left and minimum top right corners along each dimension
	RTreeNode *tnInfoMin = (RTreeNode *)malloc(DIMENSION * sizeof(RTreeNode));
	RTreeNode *tnInfoMax = (RTreeNode *)malloc(DIMENSION * sizeof(RTreeNode));
	RTreeNode temp;
	RLstNd lstNdTemp = NULL;
    int iCnt = 0;
	Boolean bIsFirst = TRUE;
	Region rgnRectTemp = NULL;
	Region rgnRectOut = NULL;
	double dNormSep;
	double dMaxNormSep;
	dimension dimMaxNormSep ;
	switch(ptrChildLst->ptrFirstNd->tnInfo->ndType)
    {	case INTNODE:   lstNdTemp = ptrChildLst->ptrFirstNd;
                        //for storing maximum bottomleft and minimum top right corners that are found till now
	                    rgnRectTemp = RinitRgnRect(NULL, NULL);
	                    // Region rgnRectTemp2= initRgnRect(NULL, NULL);
	           //          for(iCnt = 0; iCnt < DIMENSION; iCnt++)
            //             {
            //             	if(lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] > lstNdTemp->ptrNextNd->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt])
            //             	{
            //             		rgnRectTemp->iBottomLeft[iCnt] =   lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt];
            //             		rgnRectTemp->iTopRight[iCnt] = lstNdTemp->ptrNextNd->tnInfo->tdInfo->rgnRect->iTopRight[iCnt];

            //             	}
            //             	else
            //             	{
            //             		rgnRectTemp->iBottomLeft[iCnt] =   lstNdTemp->ptrNextNd->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt];
            //             		rgnRectTemp->iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt];

            //             	}
                        	 
				        // rgnRectTemp->iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt];
            //             }
	                    
                        //for storing the minim bottom left and maximum topr right corner
	                    rgnRectOut = RinitRgnRect(NULL, NULL);
                        while(lstNdTemp != NULL)
                        {     for(iCnt = 0; iCnt < DIMENSION; iCnt++)
                              {     
                              	if(bIsFirst)
                                    {     //if this is first node set its corners to max bottom left and min top right corners
				                          rgnRectTemp->iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt];
				                          rgnRectTemp->iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt];
				                          rgnRectOut->iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt];
				                          rgnRectOut->iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt];
                                          //set the marker to the at the iCnt to dimesnion to current node if it has the max bottom left and min top right corners.
				                          tnInfoMin[iCnt] = lstNdTemp->tnInfo;
				                          tnInfoMax[iCnt] = lstNdTemp->tnInfo;
				                          continue;
                                    }
			                        if(lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] > rgnRectTemp->iBottomLeft[iCnt])
                                    {     //if the current node has greater bottom left corner than the node pointed by the tnInfoMin[ iCnt ] along iCnt Dimension assign it to tnInfoMin[ iCnt ]
				                          rgnRectTemp->iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt];
				                          tnInfoMin[iCnt] = lstNdTemp->tnInfo;
                                    }
			                        if(lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt] < rgnRectTemp->iTopRight[iCnt])
                                    {     //if the current node has lower top right corner than the node pointed by the tnInfoMax[ iCnt ] along iCnt Dimension assign it to tnInfoMax[ iCnt ]
				                          rgnRectTemp->iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt];
				                          tnInfoMax[iCnt] = lstNdTemp->tnInfo;
                                    }
			                        // else
                           //          {     //make sure that the max bottom left and min top right are not for the same rectangle
				                       //    if(tnInfoMin[iCnt] == tnInfoMax[iCnt])
                           //                      tnInfoMax[ iCnt ] = lstNdTemp->tnInfo;
                           //          }
                                    //sotre the minimum bottom left corner along each dimension
			                        if(lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] < rgnRectOut->iBottomLeft[iCnt])
                                                rgnRectOut->iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt];
                                    //sotre the maximum top right corner along each dimension
			                        if(lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt] > rgnRectOut->iTopRight[iCnt])
                                                rgnRectOut->iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt];
                               }
                               lstNdTemp = lstNdTemp->ptrNextNd;
		                       if(bIsFirst)
                                    bIsFirst = FALSE;
                         }	//while
	                     dNormSep = 0;
						 dMaxNormSep = 0;
	                     dimMaxNormSep = 0;
	                     for(iCnt = 0; iCnt < DIMENSION; iCnt++)
                         {     //calculate normal seperation along each dimension
                               dNormSep = fabs(rgnRectTemp->iBottomLeft[iCnt] - rgnRectTemp->iTopRight[iCnt]) / fabs(rgnRectOut->iTopRight[iCnt] - rgnRectOut->iBottomLeft[iCnt]);
                               if(dNormSep > dMaxNormSep)
                               {   dMaxNormSep = dNormSep;
			                       dimMaxNormSep = iCnt;
                               }
                         }
                         if(tnInfoMin[(int)dimMaxNormSep] == tnInfoMax[(int)dimMaxNormSep])
                         {
                         	lstNdTemp = ptrChildLst->ptrFirstNd;
                         	temp=tnInfoMax[(int)dimMaxNormSep];
                         	if(temp != lstNdTemp->tnInfo)
                         	{
                         		rgnRectTemp->iTopRight[(int)dimMaxNormSep] = lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[(int)dimMaxNormSep];
                         		 tnInfoMax[(int)dimMaxNormSep] = lstNdTemp->tnInfo;
                         		lstNdTemp = lstNdTemp->ptrNextNd;
                         	}
                         	else
                         	{
                         		rgnRectTemp->iTopRight[(int)dimMaxNormSep] = lstNdTemp->ptrNextNd->tnInfo->tdInfo->rgnRect->iTopRight[(int)dimMaxNormSep];
                         		tnInfoMax[(int)dimMaxNormSep] = lstNdTemp->ptrNextNd->tnInfo;
                         		lstNdTemp = lstNdTemp->ptrNextNd->ptrNextNd;
                         	}
                         	while(lstNdTemp != NULL)
                         	{
                         		if(lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[(int)dimMaxNormSep] < rgnRectTemp->iTopRight[(int)dimMaxNormSep] && temp != lstNdTemp->tnInfo)
                         		{
                         			rgnRectTemp->iTopRight[(int)dimMaxNormSep] = lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[(int)dimMaxNormSep];
				                          tnInfoMax[(int)dimMaxNormSep] = lstNdTemp->tnInfo;
                         		}
                         			lstNdTemp = lstNdTemp->ptrNextNd;
                         	} 
                         }
                         if(tnInfoMin[(int)dimMaxNormSep]==tnInfoMax[(int)dimMaxNormSep])
                         {
                         	printf("error in the code\n");
                         	exit(-1);
                         }
                         // remove the node pointed by tnInfoMin at tnInfoMax at dMaxNormSep and assign them to be the two nodes that are picked for split
                         // printf("inside pickSeeds\n");
                         // printAllChildrenInList(ptrChildLst);
                         //  printf("printing first node to be deleted\n");
                         // printTnInfo(tnInfoMin[(int)dimMaxNormSep]);
	                     *lstNdChildOne = RdeleteLstElem(ptrChildLst, tnInfoMin[(int)dimMaxNormSep]);
	                     // printf("printing first node deleted\n");
	                     // printNode(*lstNdChildOne);
	                     // printf("intermediate pickSeeds\n");
                      //    printAllChildrenInList(ptrChildLst);
                      //    printf("printing second node to be deleted\n");
                      //    printTnInfo(tnInfoMax[(int)dimMaxNormSep]);
	                     // if(tnInfoMin[(int)dimMaxNormSep]==tnInfoMax[(int)dimMaxNormSep])
	                     // 		printAllChildrenInList(ptrChildLst);
	                     *lstNdChildTwo = RdeleteLstElem(ptrChildLst, tnInfoMax[(int)dimMaxNormSep]);
	                      // if(tnInfoMin[(int)dimMaxNormSep]==tnInfoMax[(int)dimMaxNormSep])
	                     	// 	printAllChildrenInList(ptrChildLst);
	     //                 printf("printing second node deleted\n");
	     //                 printNode(*lstNdChildTwo);
						// printf("after both deletion\n");
      //                    printAllChildrenInList(ptrChildLst);
	                     // if(*lstNdChildTwo==NULL)//in case both seeds pick up same child
	                     // {
	                     // 	// if(ptrChildLst->uiCnt==0)
	                     // 	// 	printf("list isLstEmpty\n");
	                     // 	// printf("lstNdTemp NULL %d\n",ptrChildLst->uiCnt);
	                     // 	// printAllChildrenInList(ptrChildLst);
	                     // 	// getchar();
	                     // 	*lstNdChildTwo = deleteLstFirst(ptrChildLst);
	                     // 	// printAllChildrenInList(ptrChildLst);
	                     // 	// getchar();
	                     // }
	                     // if(*lstNdChildTwo==NULL)
	                     // 	printf("lstNdChildTwo nul\n");
	                     free(rgnRectTemp->iBottomLeft);
	                     free(rgnRectTemp->iTopRight);
	                     free(rgnRectTemp);
	                     free(rgnRectOut->iBottomLeft);
	                     free(rgnRectOut->iTopRight);
	                     free(rgnRectOut);
                         break;

	case EXTNODE:        //same as in case of INTNODE but only with the difference that max seperation instead of normal seperation is picked.. so just find two nodes which are at the farthest distance along each dimension
	                     lstNdTemp = ptrChildLst->ptrFirstNd;
	                     DataPoint iDataMin = (DataPoint)calloc(DIMENSION, sizeof(dataPoint));
	                     DataPoint iDataMax = (DataPoint)calloc(DIMENSION, sizeof(dataPoint));
                         bIsFirst = TRUE;
	                     while(lstNdTemp != NULL)
                         {     for(iCnt = 0; iCnt < DIMENSION; iCnt++)
                               {   if(bIsFirst)
                                   {    iDataMin[iCnt] = lstNdTemp->tnInfo->tdInfo->dataClstElem->iData[iCnt];
				                        iDataMax[iCnt] = lstNdTemp->tnInfo->tdInfo->dataClstElem->iData[iCnt];
                                       	tnInfoMin[iCnt] = lstNdTemp->tnInfo;
				                        tnInfoMax[iCnt] = lstNdTemp->tnInfo;
			                         	continue;
                                   }
			                       if(lstNdTemp->tnInfo->tdInfo->dataClstElem->iData[iCnt] <= iDataMin[iCnt])
                                   {    iDataMin[iCnt] = lstNdTemp->tnInfo->tdInfo->dataClstElem->iData[iCnt];
				                        tnInfoMin[iCnt] = lstNdTemp->tnInfo;
                               	   }
		                           if(lstNdTemp->tnInfo->tdInfo->dataClstElem->iData[iCnt] > iDataMax[iCnt])
                                   {    iDataMax[iCnt] = lstNdTemp->tnInfo->tdInfo->dataClstElem->iData[iCnt];
                                        tnInfoMax[iCnt]  = lstNdTemp->tnInfo;
                                   }
                           	    }
                                lstNdTemp = lstNdTemp->ptrNextNd;
                               	if(bIsFirst)
                                    bIsFirst = FALSE;
                         } //while
	                     double dSep = 0;
	                     double dMaxSep = 0;
	                     dimension dimMaxSep = 0;
	                     for(iCnt = 0; iCnt < DIMENSION; iCnt++)
                         {      dSep = fabs(iDataMax[iCnt] - iDataMin[iCnt]);
		                        if(dSep > dMaxSep)
                                {   dMaxSep = dSep;
			                        dimMaxSep = iCnt;
	                            }
                         }
	                     *lstNdChildOne = RdeleteLstElem(ptrChildLst, tnInfoMin[(int)dimMaxSep]);
	                     *lstNdChildTwo = RdeleteLstElem(ptrChildLst, tnInfoMax[(int)dimMaxSep]);

	                     free(iDataMin);
	                     free(iDataMax);
                         break;
	}//switch

	free(tnInfoMin);
	free(tnInfoMax);

	return;
}

void RsplitNode(RLstNd ptrChild)
{    
	// splits the node into two nodes
     //printf("\nBHAIYYA SPLIT\n");
     //printf("%lf\t%lf\n",ptrChild->tnInfo->tdInfo->dataClstElem->iData[0],ptrChild->tnInfo->tdInfo->dataClstElem->iData[0]);
	if(ptrChild == NULL || RisLstEmpty(ptrChild->ptrChildLst))
		return;

	RLstNd lstNdOne = NULL;
	RLstNd lstNdTwo = NULL;
	RTreeNode tnInfoTemp = NULL;

	RLstNd lstNdTemp = NULL;

	double dExpOne = -1;
	double dExpTwo = -1;

    //pick two nodes that are farthest along any dimension
	RpickSeeds(ptrChild->ptrChildLst, &lstNdOne, &lstNdTwo);

	if(lstNdOne == NULL || lstNdTwo == NULL)
		return;
    //printf("BHAIYYA SPLIT\t");
    //printf("%lf\t%lf\n",ptrChild->tnInfo->tdInfo->dataClstElem->iData[0],ptrChild->tnInfo->tdInfo->dataClstElem->iData[1]);

    //create two child lists
	RLstNd ptrChildTemp = RinitLstNd(RinitIntNd(NULL, NULL));
	RLstNd ptrChildSib = RinitLstNd(RinitIntNd(NULL, NULL));
//	RTreeNode tn = RinitIntNd(NULL, NULL);

    //link the two child lists so that one follows the other
	ptrChildTemp->ptrChildLst = RinitHdrNd();
	ptrChildSib->ptrChildLst = RinitHdrNd();
	ptrChildSib->ptrNextNd = ptrChild->ptrNextNd;

    //insert the picked children one into each of the list
	RinsertLstNd(ptrChildTemp->ptrChildLst, lstNdOne);
	RsetRect(ptrChildTemp, lstNdOne->tnInfo);
	RinsertLstNd(ptrChildSib->ptrChildLst, lstNdTwo);
	RsetRect(ptrChildSib, lstNdTwo->tnInfo);

	Region rgnNewRectOne = RinitRgnRect(NULL, NULL);
	Region rgnNewRectTwo = RinitRgnRect(NULL, NULL);

	Boolean bIsOne = FALSE;
	Boolean bIsNdOneInComp = FALSE;
	Boolean bIsNdTwoInComp = FALSE;

	int iCnt = 0;

	lstNdTemp = RdeleteLstFirst(ptrChild->ptrChildLst);

    //pick one element from the list of children of the node to be split
	while(lstNdTemp != NULL)
    {   //if one of the nodes has so few entires that all the remaining children are to be assigned set that node to be incomplete
		if(ptrChildTemp->ptrChildLst->uiCnt + ptrChild->ptrChildLst->uiCnt == RMINENTRIES - 1)
			bIsNdOneInComp = TRUE;

		if(ptrChildSib->ptrChildLst->uiCnt + ptrChild->ptrChildLst->uiCnt == RMINENTRIES - 1)
			bIsNdTwoInComp = TRUE;
        //if both nodes are not potentiall incomplete i.e. all the remaining children need not be assigned to it for the node not to underflow
		if(!bIsNdOneInComp && !bIsNdTwoInComp)
        {   dExpOne = -1;
		    dExpTwo = -1;
            //find the area by which the rectangle in each node should be expanded to accomodat the given rectangle
		    RexpansionArea(ptrChildTemp->tnInfo->tdInfo->rgnRect, lstNdTemp->tnInfo, &dExpOne, rgnNewRectOne);
		    RexpansionArea(ptrChildSib->tnInfo->tdInfo->rgnRect, lstNdTemp->tnInfo, &dExpTwo, rgnNewRectTwo);
            //printf("kahan:%d\t%d\t%d\t%d\n",rgnNewRectOne->iBottomLeft[0], rgnNewRectOne->iBottomLeft[1], rgnNewRectOne->iTopRight[0], rgnNewRectOne->iTopRight[1]);

            //ark the node that requires least expansion to be the list into which the child can be added
		    if(dExpOne < dExpTwo)
			     bIsOne = TRUE;
	        if(dExpOne > dExpTwo)
		         bIsOne = FALSE;
	        if(dExpOne == dExpTwo)
            {    //incase both reequired to be enlarged by same amount pick the one with lower area currently
		         double dAreaOne = Rarea(ptrChildTemp->tnInfo->tdInfo->rgnRect);
			     double dAreaTwo = Rarea(ptrChildSib->tnInfo->tdInfo->rgnRect);
			     if(dAreaOne < dAreaTwo)
                      bIsOne = TRUE;
                 if(dAreaOne > dAreaTwo)
                      bIsOne = FALSE;
                 if(dAreaOne == dAreaTwo)
                 {    //incase the area are same too pick the node which has lesser number of children
                      if(ptrChildTemp->ptrChildLst->uiCnt < ptrChildSib->ptrChildLst->uiCnt)
                           bIsOne = TRUE;
                      else
					       bIsOne = FALSE;
                 }
           }
		} //if
		
		else
        {   //if one of the nodes is potentially incomplete mark it to be the node to which the child has to be assigned
		    if(bIsNdOneInComp)
                  bIsOne = TRUE;
		    if(bIsNdTwoInComp)
			      bIsOne = FALSE;
		}

		if(bIsOne)
        {   //insert the new child
			RinsertLstNd(ptrChildTemp->ptrChildLst, lstNdTemp);
            if(bIsNdOneInComp)
            {   dExpOne = -1;
				RexpansionArea(ptrChildTemp->tnInfo->tdInfo->rgnRect, lstNdTemp->tnInfo, &dExpOne, rgnNewRectOne);
			}
            //expand the rectangle to accomodate new child
			for(iCnt = 0; iCnt < DIMENSION; iCnt++)
            {	ptrChildTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] = rgnNewRectOne->iBottomLeft[iCnt];
				ptrChildTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt] = rgnNewRectOne->iTopRight[iCnt];
			}
		}
		else
        {   //insert the new child
			RinsertLstNd(ptrChildSib->ptrChildLst, lstNdTemp);
            if(bIsNdTwoInComp)
            {   dExpTwo = -1;
				RexpansionArea(ptrChildSib->tnInfo->tdInfo->rgnRect, lstNdTemp->tnInfo, &dExpTwo, rgnNewRectTwo);
			}
            //expand the rectangle to accomodate the new child
			for(iCnt = 0; iCnt < DIMENSION; iCnt++)
            {	ptrChildSib->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] = rgnNewRectTwo->iBottomLeft[iCnt];
				ptrChildSib->tnInfo->tdInfo->rgnRect->iTopRight[iCnt] = rgnNewRectTwo->iTopRight[iCnt];
			}
		}
        //picke next node
		lstNdTemp = RdeleteLstFirst(ptrChild->ptrChildLst);
	}
    //set the node that is passed too first of the two nodes and set the next pointer of the second to the next pointer of the node that is passed..so that now in place of the node that is passed we have to nodes instead
	ptrChildTemp->ptrChildLst->ptrParentNd = ptrChildTemp;
	ptrChildSib->ptrChildLst->ptrParentNd = ptrChildSib;
	ptrChildTemp->ptrNextNd = ptrChildSib;
	ptrChild->ptrChildLst->ptrParentNd = NULL;
	free(ptrChild->ptrChildLst);

	free(ptrChild->tnInfo->tdInfo->rgnRect->iBottomLeft);
	free(ptrChild->tnInfo->tdInfo->rgnRect->iTopRight);
	free(ptrChild->tnInfo->tdInfo->rgnRect);

	free(ptrChild->tnInfo->tdInfo);
	free(ptrChild->tnInfo);

	ptrChild->tnInfo = ptrChildTemp->tnInfo;
	ptrChild->ptrChildLst = ptrChildTemp->ptrChildLst;
	ptrChild->ptrNextNd = ptrChildTemp->ptrNextNd;

	ptrChildTemp->ptrNextNd = NULL;
	ptrChildTemp->ptrChildLst = NULL;
	ptrChildTemp->tnInfo = NULL;
	free(ptrChildTemp);

	free(rgnNewRectOne->iBottomLeft);
	free(rgnNewRectOne->iTopRight);
	free(rgnNewRectOne);
	free(rgnNewRectTwo->iBottomLeft);
	free(rgnNewRectTwo->iTopRight);
	free(rgnNewRectTwo);

	return;
}


Boolean RinsertTree(RHdrNd hdrNdTree, RTreeNode tnInfo)
{   //insert a node into tree
	int iCnt = 0;

	if(hdrNdTree->ptrFirstNd == NULL || hdrNdTree->ptrFirstNd->tnInfo->ndType == EXTNODE)
    {   
     //incase of the external node insert the node into the child list and if the node has more tha max entries return true to indicat that
		// if(tnInfo->tdInfo->dataClstElem->iNum ==102 )
		// 	{
		// 		printf("inserting %d %lf %lf finally at the leaf node\n",tnInfo->tdInfo->dataClstElem->iNum,tnInfo->tdInfo->dataClstElem->iData[0],tnInfo->tdInfo->dataClstElem->iData[1] );	
		// 		RprintTree(hdrNdTree);
		// 		getchar();
		// 	}
		
		RinsertLstElem(hdrNdTree, tnInfo);
		//should be there: printf(".... ...entries nw after insertion is : %d and limit is %d\n",hdrNdTree->uiCnt,RMAXENTRIES);
		// if(tnInfo->tdInfo->dataClstElem->iNum ==102)
		// 	{
		// 		printf("after inserting %d %lf %lf finally at the leaf node\n",tnInfo->tdInfo->dataClstElem->iNum,tnInfo->tdInfo->dataClstElem->iData[0],tnInfo->tdInfo->dataClstElem->iData[1] );	
		// 		RprintTree(hdrNdTree);
		// 		printf("size of hdrNdTree %d max %d\n",hdrNdTree->uiCnt,RMAXENTRIES );
		// 		getchar();
		// 	}
		
		return(hdrNdTree->uiCnt > RMAXENTRIES) ? TRUE : FALSE;
    }

	// if(tnInfo->tdInfo->dataClstElem->iNum ==102)
	// {
	// 	printf("before setting rectangle %d %lf %lf finally at the leaf node\n",tnInfo->tdInfo->dataClstElem->iNum,tnInfo->tdInfo->dataClstElem->iData[0],tnInfo->tdInfo->dataClstElem->iData[1] );	
	// 	getchar();
	// 	RprintTree(hdrNdTree);
	// 	getchar();
	// }
	if(hdrNdTree->ptrFirstNd->ptrChildLst->uiCnt == 0)
		RsetRect(hdrNdTree->ptrFirstNd, tnInfo);

// if(tnInfo->tdInfo->dataClstElem->iNum ==102)
// 	{
// 		printf("after setting rectangle and before pick child %d %lf %lf finally at the leaf node\n",tnInfo->tdInfo->dataClstElem->iNum,tnInfo->tdInfo->dataClstElem->iData[0],tnInfo->tdInfo->dataClstElem->iData[1] );	
// 		getchar();
// 		RprintTree(hdrNdTree);
// 		getchar();
// 	}

    //pick the child into which the node can be inserted
   
	RLstNd lstNdTemp = RpickChild(hdrNdTree, tnInfo);
	//expnandRect( lstNdTemp, tnInfo );
	// if(tnInfo->tdInfo->dataClstElem->iNum ==102)
	// {
	// 	printf("after pick child and before insert %d %lf %lf finally at the leaf node\n",tnInfo->tdInfo->dataClstElem->iNum,tnInfo->tdInfo->dataClstElem->iData[0],tnInfo->tdInfo->dataClstElem->iData[1] );	
	// 	getchar();
	// 	RprintTree(hdrNdTree);
	// 	printf("\n\nafter pick child %lf %lf %lf %lf\n",lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[0],lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[1],lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[0],lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[1] );
	// 	getchar();
	// }
	

   //call insert tree on the child that is picked
	if(RinsertTree(lstNdTemp->ptrChildLst, tnInfo))
    {   //incase the resultant insert has caused the node to over flow invoke split node
        //should be there: printf("\nBHAIYYA SPLIT\n");
        //should be there: printf("%lf\t%lf\t%d\n",lstNdTemp->tnInfo->tdInfo->dataClstElem->iData[0],lstNdTemp->tnInfo->tdInfo->dataClstElem->iData[0],isLstEmpty(lstNdTemp->ptrChildLst));
  //       if(tnInfo->tdInfo->dataClstElem->iNum ==102)
		// {
		// 	printf("before split node %d %lf %lf\n",tnInfo->tdInfo->dataClstElem->iNum,tnInfo->tdInfo->dataClstElem->iData[0],tnInfo->tdInfo->dataClstElem->iData[1] );	
		// 	getchar();
		// 	RprintTree(hdrNdTree);
		// 	getchar();
		// }
        RsplitNode(lstNdTemp);
  //       if(tnInfo->tdInfo->dataClstElem->iNum ==102)
		// {
		// 	printf("after split node %d %lf %lf\n",tnInfo->tdInfo->dataClstElem->iNum,tnInfo->tdInfo->dataClstElem->iData[0],tnInfo->tdInfo->dataClstElem->iData[1] );	
		// 	getchar();
		// 	RprintTree(hdrNdTree);
		// 	getchar();
		// }
	    hdrNdTree->uiCnt++;
	    //getchar();
        //if after split node is invoked the node is overflowing return treu to to its parent to let it know that node is over flowing
	    return (hdrNdTree->uiCnt > RMAXENTRIES) ? TRUE : FALSE;
	}
	// if(tnInfo->tdInfo->dataClstElem->iNum ==102)
	// 	{
	// 		printf("end of function %d %lf %lf finally at the leaf node\n",tnInfo->tdInfo->dataClstElem->iNum,tnInfo->tdInfo->dataClstElem->iData[0],tnInfo->tdInfo->dataClstElem->iData[1] );	
	// 		getchar();
	// 		RprintTree(hdrNdTree);
	// 		getchar();
	// 	}

	return FALSE;
}

RHdrNd RbuildRTree(DataHdr dataHdrLst)
{   // build RTree with the list of elements passed
	RHdrNd hdrNdTree = RinitHdrNd();
	hdrNdTree->ptrFirstNd = RinitLstNd(RinitIntNd(NULL, NULL));
	hdrNdTree->uiCnt++;
	hdrNdTree->ptrFirstNd->ptrChildLst = RinitHdrNd();

	//DataNd dataTemp = dataHdrLst->dataFirstNd;
    ////should be there: printf("MIn %d max %d\n",MINENTRIES,RMAXENTRIES);
	int cnt = 0,i;
    //retirve element one by one and insert them into tree invoking create root incase of root split
	for(i=0;i<dataHdrLst->uiCnt;i++)
	{
		//printf("i=%d\n",i );
		RinsertTree(hdrNdTree, RinitExtNd((dataHdrLst->dataClstElem)+i));
		if(hdrNdTree->uiCnt > 1)
			hdrNdTree = RcreateRoot(hdrNdTree);
	}

	return	hdrNdTree;
}

RHdrNd RbuildRTreeFromCells(BCellListHd cellsList, Region epsExtendedRgn)
{
	RHdrNd hdrNdTree = RinitHdrNd();
	hdrNdTree->ptrFirstNd = RinitLstNd(RinitIntNd(NULL, NULL));
	hdrNdTree->uiCnt++;
	hdrNdTree->ptrFirstNd->ptrChildLst = RinitHdrNd();

	//DataNd dataTemp = dataHdrLst->dataFirstNd;
    ////should be there: printf("MIn %d max %d\n",RMINENTRIES,RMAXENTRIES);
	int cnt = 0,i;
	int j=0;
	Boolean cellFlag;
    //retrieve element one by one and insert them into tree invoking create root incase of root split
	
	BCellListNode currCellNode = cellsList->first;

	for(i=0;i<cellsList->count;i++)
	{
		// BCell currCell = currCellNode->bCellElem;
		// cellFlag = GisContainsCell(epsExtendedRgn,currCell);

		// if(cellFlag==TRUE)
		// {
		// 	CellDataHd currCellDataList = currCell->cellDataList;
		// 	CellData currCellData = currCellDataList->first;
			
		// 	while(currCellData!=NULL)
		// 	{
		// 		/*printf("inserting %d %lf %lf\n",currCellData->data->iNum, currCellData->data->iData[0], currCellData->data->iData[1] );
		// 		if(currCellData->data->iNum==58 || currCellData->data->iNum==62 || currCellData->data->iNum==63)
		// 			getchar(); */
				
		// 		RinsertTree(hdrNdTree, RinitExtNd(currCellData->data));
		// 		//printf("root count = %d\n",hdrNdTree->uiCnt );
		// 		// if(currCellData->data->iNum==58 || currCellData->data->iNum==62 || currCellData->data->iNum==63)
		// 		// {
		// 		// 	getchar();
		// 		// 	RprintTree(hdrNdTree);
		// 		// }
					
		// 		if(hdrNdTree->uiCnt > 1)
		// 			hdrNdTree = RcreateRoot(hdrNdTree);

		// 		// if(currCellData->data->iNum==58 || currCellData->data->iNum==62 || currCellData->data->iNum==63)
		// 		// {
		// 		// 	getchar();
		// 		// 	RprintTree(hdrNdTree);
		// 		// }

		// 		currCellData=currCellData->next;			

		// 	}
		// }
		// else
		// {
		// 	CellDataHd currCellDataList = currCell->cellDataList;
		// 	CellData currCellData = currCellDataList->first;
			
		// 	while(currCellData!=NULL)
		// 	{
		// 		/*printf("inserting %d %lf %lf\n",currCellData->data->iNum, currCellData->data->iData[0], currCellData->data->iData[1] );
		// 		if(currCellData->data->iNum==58 || currCellData->data->iNum==62 || currCellData->data->iNum==63)
		// 			getchar(); */

		// 		Data currData = currCellData->data;
		// 		if(RisContains(epsExtendedRgn,currData))
		// 		{
		// 			RinsertTree(hdrNdTree, RinitExtNd(currCellData->data));
		// 			//printf("root count = %d\n",hdrNdTree->uiCnt );
		// 			// if(currCellData->data->iNum==58 || currCellData->data->iNum==62 || currCellData->data->iNum==63)
		// 			// {
		// 			// 	getchar();
		// 			// 	RprintTree(hdrNdTree);
		// 			// }
						
		// 			if(hdrNdTree->uiCnt > 1)
		// 				hdrNdTree = RcreateRoot(hdrNdTree);

		// 			// if(currCellData->data->iNum==58 || currCellData->data->iNum==62 || currCellData->data->iNum==63)
		// 			// {
		// 			// 	getchar();
		// 			// 	RprintTree(hdrNdTree);
		// 			// }
		// 		}
		// 		currCellData=currCellData->next;					

		// 	}	
		
		// }


			BCell currCell = currCellNode->bCellElem;			
			CellDataHd currCellDataList = currCell->cellDataList;
			CellData currCellData = currCellDataList->first;
			
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

		currCellNode = currCellNode->next;
		
	}

	return	hdrNdTree;
}

void RprintTree(RHdrNd hdrNdTree)
{   //travers along the tree and print the tree
	RLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
	int iCnt = 0;
	static int iIndent = 0;
	iIndent++;

	while(lstNdTemp != NULL)
    {   for(iCnt = 0; iCnt < iIndent; iCnt++)
			printf("---");

		if(lstNdTemp->tnInfo->ndType == INTNODE)
        {   printf("i hav %d children..\n",lstNdTemp->ptrChildLst->uiCnt);
			printf(" BottomLeft: ");
			for(iCnt = 0; iCnt < DIMENSION; iCnt++)
				printf("%lf ", lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt]);
            printf(" TopRight: ");
			for(iCnt = 0; iCnt < DIMENSION; iCnt++)
				printf("%lf ", lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt]);
            printf("\n");
			printf("i hav %d children..\n",lstNdTemp->ptrChildLst->uiCnt);
			RprintTree(lstNdTemp->ptrChildLst);
		}
		else
        {   printf(" DataPoint: ");
    		printData(lstNdTemp->tnInfo->tdInfo->dataClstElem);
			// for(iCnt = 0; iCnt < DIMENSION; iCnt++)
			// 	printf("%lf ", lstNdTemp->tnInfo->tdInfo->dataClstElem->iData[iCnt]);
			printf("\n");
		}
		lstNdTemp = lstNdTemp->ptrNextNd;
	}
	iIndent--;
	return;
}

double RfindDist(DataPoint iDataOne, DataPoint iDataTwo)
{   
	double dDist = 0;
	int iCnt = 0;
	int k=0;
	for(iCnt = 0; iCnt < DIMENSION; iCnt++)
	{    
		dDist = dDist + pow(iDataOne[iCnt] - iDataTwo[iCnt], 2);
	}

	return sqrt(dDist);
}

unsigned int RgetNeighborHood(RHdrNd hdrNdTree, Data dataTemp,int size_data)
{   
	
	if(hdrNdTree == NULL || dataTemp == NULL)
		return 0;

	Region rgnRect = RinitRgnRect(NULL, NULL);
	int iCnt = 0;
	unsigned int uiRecCnt = 0;
	for(iCnt = 0; iCnt < DIMENSION; iCnt++)
    {   
		rgnRect->iBottomLeft[iCnt] = dataTemp->iData[iCnt] - EPS;
		rgnRect->iTopRight[iCnt] = dataTemp->iData[iCnt] + EPS;
	}
	dataTemp->neighbors = RinitNbHdr();
	
	assert(dataTemp->neighbors!=NULL);
	//printf("%d No of neighbors %d\n", dataTemp->iNum, dataTemp->neighbors->nbhCnt);
	//printf("\nPrinting the datapoint for which neighborhood is being queried\n");
	//printData(dataTemp);	
	//printf("\nPrinting the eps extended region for which neighborhood is being queried\n");
	//GprintRegion(rgnRect);

    uiRecCnt = RfindRecords(hdrNdTree, rgnRect, dataTemp);

    // if(dataTemp->neighbors->nbhCnt >= MINPOINTS)
    // {   
    // 	dataTemp->core_tag = TRUE;
    // } 
	free(rgnRect->iBottomLeft);
	free(rgnRect->iTopRight);
	free(rgnRect);

	return uiRecCnt;
}

unsigned int RfindRecords(RHdrNd hdrNdTree, Region rgnRect, Data dataTemp)
{   
	//finds the record in the given search rectangle and returns  number of records found
    //if datapoint is passed in iData finds the number of records in eps neighborhood

	if(RisLstEmpty(hdrNdTree) || rgnRect == NULL || rgnRect->iBottomLeft == NULL || rgnRect->iTopRight == NULL)
		return 0;

	unsigned int uiRecCnt = 0;
	int iCnt = 0;
    double t;
    static int flag = 0;
	RLstNd lstNdTemp = hdrNdTree->ptrFirstNd;

	while(lstNdTemp != NULL)
    {	//nodesVisited++;
		switch(lstNdTemp->tnInfo->ndType)
        {   case INTNODE:   //incase of internal node if the node is over lapping with search rectangle descend into the node and add to the count the number of new records found if any
				            if(RisOverLap(lstNdTemp->tnInfo->tdInfo->rgnRect, rgnRect))
                                   uiRecCnt += RfindRecords(lstNdTemp->ptrChildLst, rgnRect, dataTemp);
                            //printf("\nuiRecCnt = %d", uiRecCnt);
                            break;
            case EXTNODE:   //incase of external node if the child node is contanied in serach rectangle increment the count and return the count to the parent
				            if(RisContains(rgnRect, lstNdTemp->tnInfo->tdInfo->dataClstElem->iData))
                            {      if(dataTemp->iData == NULL)
                                   {     //for( iCnt = 0; iCnt < DIMENSION; iCnt++ )
						                 //	printf( "%d ", lstNdTemp->tnInfo->tdInfo->dataClstElem->iData[ iCnt ] );
			                             //printf( "\n" );
                                   		//printf("%d\n",lstNdTemp->tnInfo->tdInfo->dataClstElem->iNum );
                                         RinsertNeighbors(dataTemp, lstNdTemp->tnInfo->tdInfo->dataClstElem, 0);
						                 uiRecCnt++;
                                   }
					               else
                                   {      //see if the node is with in the esp neighborhood and increment the count if it is
                                          if(lstNdTemp->tnInfo->tdInfo->dataClstElem->iData != dataTemp->iData)
                                          {     //double t;
							                    if((t = RfindDist(lstNdTemp->tnInfo->tdInfo->dataClstElem->iData, dataTemp->iData)) <= EPS)
                                                {
		                                        /*for( iCnt = 0; iCnt < DIMENSION; iCnt++ )
							                    printf("nbg %.0f ", lstNdTemp->tnInfo->tdInfo->dataClstElem->iData[iCnt]);
						                        printf("\n");*/
						                        //printf("\nDistance = %lf", t);
								                //insertDataLstNd(dataSeeds, lstNdTemp->tnInfo->tdInfo->dataClstElem);
								                //printf("%d\n",lstNdTemp->tnInfo->tdInfo->dataClstElem->iNum );
								                RinsertNeighbors(dataTemp, lstNdTemp->tnInfo->tdInfo->dataClstElem, t);
								                uiRecCnt++;
								                //printf("\nuiRecCnt = %d", uiRecCnt);
								                //printf("\nuiRecCnt : %d", uiRecCnt);
			                                    }
                                          }
		                           }
                             }
	                         break;

		}//switch
		lstNdTemp = lstNdTemp->ptrNextNd;
	}  
	return uiRecCnt;
}

Boolean RisContains(Region rgnRect, DataPoint iData)
{   //return TRUE if the data point is contained in a rectangle
	int iCnt = 0;
	Boolean bIsContains = TRUE;
	for(iCnt = 0; iCnt < DIMENSION; iCnt++)
    {   if((rgnRect->iBottomLeft[iCnt] > iData[iCnt]) || (rgnRect->iTopRight[iCnt] < iData[iCnt]))
        {   bIsContains = FALSE;
			break;
		}
	}
	return bIsContains;
}

Boolean RisOverLap(Region rgnRectOne, Region rgnRectTwo)
{   //returs TRUE if the two rectangles are over lapping
	int iDim = 0;
	for(iDim = 0; iDim < DIMENSION; iDim++)
		if(rgnRectOne->iTopRight[iDim] < rgnRectTwo->iBottomLeft[iDim] || rgnRectTwo->iTopRight[iDim] < rgnRectOne->iBottomLeft[iDim])
			return FALSE;
	return TRUE;
}


void RappendRTree(RHdrNd hdrNdTree, DataHdr dataHdrLst)
{   // build RTree with the list of elements passed
	
    //DataNd dataTemp = dataHdrLst->dataFirstNd;
    //should be there: printf("MIn %d max %d\n",RMINENTRIES,RMAXENTRIES);
	int cnt = 0,i;
    //retirve element one by one and insert them into tree invoking create root incase of root split
    for(i=0;i<dataHdrLst->uiCnt;i++)
	{
		RinsertTree(hdrNdTree, RinitExtNd(dataHdrLst->dataClstElem+i));
	}

	return;	
}


void freeRTree(RHdrNd hdrNdTree)
{
	// code to Free R Tree	

	if(hdrNdTree == NULL) 
		return;

	if(hdrNdTree->uiCnt == 0)
	{
		free(hdrNdTree);
		return;
	}
	
	RLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
	RLstNd lstNdNextTemp;
	
	if(lstNdTemp!=NULL)
	{
		while(lstNdTemp != NULL)
    	{
			switch(lstNdTemp->tnInfo->ndType)
        	{   case INTNODE:   
				    			free(lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft);
				    			free(lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight);
				    			free(lstNdTemp->tnInfo->tdInfo->rgnRect);
				    			free(lstNdTemp->tnInfo->tdInfo);
				    			free(lstNdTemp->tnInfo);				    			
				    			freeRTree(lstNdTemp->ptrChildLst);                            	
                            	lstNdNextTemp = lstNdTemp->ptrNextNd;
                            	free(lstNdTemp);
                            	break;
            	case EXTNODE:   					            
					            free(lstNdTemp->tnInfo->tdInfo);
					            free(lstNdTemp->tnInfo);
					            freeRTree(lstNdTemp->ptrChildLst);
					            lstNdNextTemp = lstNdTemp->ptrNextNd;
                            	free(lstNdTemp);	                            
		                        break;

			}//switch
			
			lstNdTemp = lstNdNextTemp;
		}
	}	
	
	free(hdrNdTree);	
	return;

}

void isCorrectRTree(RHdrNd hdrNdTree)
{
    RLstNd lstNdTemp = hdrNdTree->ptrFirstNd,temp;
    int iCnt = 0,flag;

    while(lstNdTemp != NULL)
        {   
        if(lstNdTemp->tnInfo->ndType == INTNODE)
            {   
            temp=lstNdTemp->ptrChildLst->ptrFirstNd;
            while(temp!=NULL)
            {
                flag=0;
                for(iCnt=0;iCnt<DIMENSION;iCnt++)
                {
                    if(temp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt]<lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt])
                        flag=1;
                }
                for(iCnt=0;iCnt<DIMENSION;iCnt++)
                {
                    if(temp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt]>lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt])
                        flag=1;
                }
                if(flag==1)
                {
                    printf("WRONG!!!\n");
                    for(iCnt=0;iCnt<DIMENSION;iCnt++)
                    {
                        printf("%lf %lf %lf %lf\n",temp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt],temp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt],lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt],lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt]);
                    }
                }
                temp=temp->ptrNextNd;
            }
            isCorrectRTree(lstNdTemp->ptrChildLst);
        }
        lstNdTemp = lstNdTemp->ptrNextNd;
    }
    return;
}	
