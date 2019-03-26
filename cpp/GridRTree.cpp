#include "GridRTree.h"

#include <stdbool.h>
#include <vector>
#include <iostream>
using namespace std;

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

Boolean GisOverLap(Region rgnRectOne, Region rgnRectTwo)
{   //returs TRUE if the two rectangles are over lapping
	int iDim = 0;
	for(iDim = 0; iDim < DIMENSION; iDim++)
		if(rgnRectOne->iTopRight[iDim] <= rgnRectTwo->iBottomLeft[iDim] || rgnRectTwo->iTopRight[iDim] < rgnRectOne->iBottomLeft[iDim])
			return FALSE;
	return TRUE;
}

Boolean GisOverLapTotal(Region rgnRectOne, Region rgnRectTwo)
{
	// returns TRUE ONLY in case of TOTAL overlap 
	int iDim=0;
	for(iDim=0; iDim < DIMENSION; iDim++)
	{
		if(rgnRectOne->iBottomLeft[iDim] < rgnRectTwo->iBottomLeft[iDim] || rgnRectOne->iTopRight[iDim] > rgnRectTwo->iTopRight[iDim]  )
			return FALSE;
	}
	return TRUE;
}

Boolean GisOverLap2(Region rgnRectOne, Region rgnRectTwo)
{   //returs TRUE if the two rectangles are over lapping
	int iDim = 0;
	for(iDim = 0; iDim < DIMENSION; iDim++)
		//if(rgnRectOne->iTopRight[iDim] <= rgnRectTwo->iBottomLeft[iDim] || rgnRectTwo->iTopRight[iDim] < rgnRectOne->iBottomLeft[iDim])
		if(rgnRectOne->iTopRight[iDim] <= rgnRectTwo->iBottomLeft[iDim] || rgnRectTwo->iTopRight[iDim] <= rgnRectOne->iBottomLeft[iDim])
			return FALSE;
	return TRUE;
}

Region GinitRgnRect(Dimension iBottomLeft, Dimension iTopRight)
{   //initializes the rectangle with the given bottom left and top right corners
    //if the values for the corners are specified NULL, initializes a rectangle with origin as co-ordinates for both corners.
	Region rgnRect = (Region)malloc(sizeof(*rgnRect));
	assert(rgnRect!=NULL);
    
    if(iBottomLeft != NULL)
		rgnRect->iBottomLeft = iBottomLeft;
	else{
		rgnRect->iBottomLeft = (Dimension)malloc(sizeof(double)*DIMENSION);
		assert(rgnRect->iBottomLeft!=NULL);
	}
    //rgnRect->iBottomLeft = (Dimension) calloc( DIMENSION, sizeof( dimension ) );


	if(iTopRight != NULL)
		rgnRect->iTopRight = iTopRight;
	else{
		rgnRect->iTopRight = (Dimension)malloc(sizeof(double)*DIMENSION);
		assert(rgnRect->iTopRight!=NULL);

	}

	return rgnRect;
}

GTreeNode GinitIntNd(Dimension iBottomLeft, Dimension iTopRight)
{   //intializes internal node of a Tree with rectangle whose bottom left and topright corners are given

	Region rgnRect = GinitRgnRect(iBottomLeft, iTopRight);
	// initializes a rectangle with the given coordonates for the bottom left and top right corners

	GTreeData tdInfo = (GTreeData)malloc(sizeof(*tdInfo));
	assert(tdInfo!=NULL);

	
	tdInfo->rgnRect = rgnRect;
    GTreeNode tnInfo = (GTreeNode)malloc(sizeof(*tnInfo));
    assert(tnInfo!=NULL);

	
	tnInfo->ndType = INTNODE;
	tnInfo->tdInfo = tdInfo;

	return tnInfo;
}

GTreeNode GinitExtNd(Group groupElem)
{
	assert(groupElem!=NULL);

	GTreeNode tnInfo = (GTreeNode)malloc(sizeof(*tnInfo));
	assert(tnInfo!=NULL);
	
	GTreeData tdInfo = (GTreeData)malloc(sizeof(*tdInfo));
    assert(tdInfo!=NULL);

	tdInfo->groupElem = groupElem;	// Data
    tnInfo->ndType = EXTNODE;	// external node
	tnInfo->tdInfo = tdInfo;

	return tnInfo;
}

void GsetRect(GLstNd lstNd, GTreeNode tnInfo)
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
        {   
        	lstNd->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] = tnInfo->tdInfo->groupElem->bcell->minOriginalBoundary[iCnt];
			lstNd->tnInfo->tdInfo->rgnRect->iTopRight[iCnt] = tnInfo->tdInfo->groupElem->bcell->maxOriginalBoundary[iCnt];
		}
		break;
	}

	return;
}

GLstNd GpickChild(GHdrNd ptrChildLst, GTreeNode tnInfo)
{  // decides which node among the child nodes to be picked for insertion and returns a pointer to that node
    if(ptrChildLst == NULL)
		return NULL;

	GLstNd lstNdTemp = ptrChildLst->ptrFirstNd;
	GLstNd lstNdChild = NULL;
	double dMinExp = -1;
    Region rgnNewRect = GinitRgnRect(NULL, NULL);
    Region rgnFinalRect = GinitRgnRect(NULL, NULL);
    int iCnt;

    // for each child child in the list of child nodes do the following
	while(lstNdTemp != NULL)
    {   //call the expansionArea function to determine the are by which the child node has to enlarged to accomodate the new point or region.
		if(GexpansionArea(lstNdTemp->tnInfo->tdInfo->rgnRect, tnInfo, &dMinExp, rgnNewRect))
        {//if the expansionArea return true mark the node to be the one that might be picked. if the expansion is same as one of the previous nodes then compare the ares of the current noe and the previous node and pick the one with least area.
               if(dMinExp < 0)
               {     dMinExp = 0 - dMinExp;
                     if(Garea(lstNdChild->tnInfo->tdInfo->rgnRect) > Garea(lstNdTemp->tnInfo->tdInfo->rgnRect))
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
		}
    	lstNdTemp = lstNdTemp->ptrNextNd;
	}
    //for the node that is picked assign the region pointed by new rectangle region and return the node
    Region rgnRectTemp = lstNdChild->tnInfo->tdInfo->rgnRect;
	lstNdChild->tnInfo->tdInfo->rgnRect = rgnFinalRect;

	free(rgnRectTemp->iBottomLeft);
	free(rgnRectTemp->iTopRight);
	free(rgnRectTemp);

	free(rgnNewRect->iBottomLeft);
	free(rgnNewRect->iTopRight);
	free(rgnNewRect);

	return lstNdChild;
}


Boolean GexpansionArea(Region rgnRect, GTreeNode tnInfo, Double ptrDMinExp, Region rgnNewRect)
{   // calculates if the area by which the rgnRect should be enlarged so as to include the tnInfo is less than the value pointed by ptrDMinExp and return TRUE and assigns rgnNewRect with the new enlarged rectangle.
    int iCnt = 0;
    Region rgnRectTemp = GinitRgnRect(NULL, NULL);
    for(iCnt = 0; iCnt < DIMENSION; iCnt++)
    {   switch(tnInfo->ndType)
        {   case INTNODE:
            //assign least of bottom left corner along each dimension to rgnRectTemp
			rgnRectTemp->iTopRight[iCnt] = (tnInfo->tdInfo->rgnRect->iTopRight[iCnt] > rgnRect->iTopRight[iCnt]) ? tnInfo->tdInfo->rgnRect->iTopRight[iCnt] : rgnRect->iTopRight[iCnt];
			rgnRectTemp->iBottomLeft[iCnt] = (tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] < rgnRect->iBottomLeft[iCnt]) ? tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] : rgnRect->iBottomLeft[iCnt];
            break;

		    case EXTNODE:
           //assign maximum of top right corner along each dimension to rgnRectTemp
			rgnRectTemp->iTopRight[iCnt] = (tnInfo->tdInfo->groupElem->bcell->maxOriginalBoundary[iCnt] > rgnRect->iTopRight[iCnt]) ? tnInfo->tdInfo->groupElem->bcell->maxOriginalBoundary[iCnt] : rgnRect->iTopRight[iCnt];
			rgnRectTemp->iBottomLeft[iCnt] = (tnInfo->tdInfo->groupElem->bcell->minOriginalBoundary[iCnt] < rgnRect->iBottomLeft[iCnt]) ? tnInfo->tdInfo->groupElem->bcell->minOriginalBoundary[iCnt] : rgnRect->iBottomLeft[iCnt];
			break;
		}
	}
    //calculate the difference in area for new rectangle and old rectangle
	double dExp = Garea(rgnRectTemp) - Garea(rgnRect);
    //in case there no value to compare ( -1 ) or incase the value is less than the value to be comparedcopy the rgnRectTemp to rgnRectNew to Return it.
	if(*ptrDMinExp == -1 || dExp <= *ptrDMinExp)
    {   if(dExp == *ptrDMinExp)
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

double Garea(Region rgnRect)
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

Boolean GinsertTree(GHdrNd hdrNdTree, GTreeNode tnInfo, int gMinEntries, int gMaxEntries)
{   //insert a node into tree
	int iCnt = 0;

	if(hdrNdTree->ptrFirstNd == NULL || hdrNdTree->ptrFirstNd->tnInfo->ndType == EXTNODE)
    {    //incase of the external node insert the node into the child list and if the node has more tha max entries return true to indicat that
		GinsertLstElem(hdrNdTree, tnInfo);
		//should be there: printf(".... ...entries nw after insertion is : %d and limit is %d\n",hdrNdTree->uiCnt,GMAXENTRIES);
		return(hdrNdTree->uiCnt > gMaxEntries) ? TRUE : FALSE;
    }

	if(hdrNdTree->ptrFirstNd->ptrChildLst->uiCnt == 0)
		GsetRect(hdrNdTree->ptrFirstNd, tnInfo);

    //pick the child into which the node can be inserted
	GLstNd lstNdTemp = GpickChild(hdrNdTree, tnInfo);

   //call insert tree on the child that is picked
	if(GinsertTree(lstNdTemp->ptrChildLst, tnInfo, gMinEntries, gMaxEntries))
    {   //incase the resultant insert has caused the node to over flow invoke split node
        //should be there: printf("\nBHAIYYA SPLIT\n");
        //should be there: printf("%lf\t%lf\t%d\n",lstNdTemp->tnInfo->tdInfo->dataClstElem->iData[0],lstNdTemp->tnInfo->tdInfo->dataClstElem->iData[0],isLstEmpty(lstNdTemp->ptrChildLst));
        GsplitNode(lstNdTemp,gMinEntries);
	    hdrNdTree->uiCnt++;
        //if after split node is invoked the node is overflowing return treu to to its parent to let it know that node is over flowing
	    return (hdrNdTree->uiCnt > gMaxEntries) ? TRUE : FALSE;
	}

	return FALSE;
}

GHdrNd GcreateRoot(GHdrNd hdrNdTree)
{  
	//in case of root split this is called to create a new root
    GHdrNd hdrNdRoot = GinitHdrNd();
    Dimension iBottomLeft = (Dimension)calloc(DIMENSION, sizeof(dimension));
	Dimension iTopRight = (Dimension)calloc(DIMENSION,sizeof(dimension));

	GLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
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
	hdrNdRoot->ptrFirstNd = GinitLstNd(GinitIntNd(iBottomLeft, iTopRight));
    hdrNdRoot->ptrFirstNd->ptrChildLst = hdrNdTree;
	hdrNdRoot->uiCnt = 1;

	return hdrNdRoot;
}

void GsplitNode(GLstNd ptrChild, int gminEntries)
{    
	// splits the node into two nodes
	if(ptrChild == NULL || GisLstEmpty(ptrChild->ptrChildLst))
		return;

	GLstNd lstNdOne = NULL;
	GLstNd lstNdTwo = NULL;
	GTreeNode tnInfoTemp = NULL;

	GLstNd lstNdTemp = NULL;

	double dExpOne = -1;
	double dExpTwo = -1;

    //pick two nodes that are farthest along any dimension
	GpickSeeds(ptrChild->ptrChildLst, &lstNdOne, &lstNdTwo);

	if(lstNdOne == NULL || lstNdTwo == NULL)
		return;

    //create two child lists
	GLstNd ptrChildTemp = GinitLstNd(GinitIntNd(NULL, NULL));
	GLstNd ptrChildSib = GinitLstNd(GinitIntNd(NULL, NULL));
	//GTreeNode tn = GinitIntNd(NULL, NULL);

    //link the two child lists so that one follows the other
	ptrChildTemp->ptrChildLst = GinitHdrNd();
	ptrChildSib->ptrChildLst = GinitHdrNd();
	ptrChildSib->ptrNextNd = ptrChild->ptrNextNd;

    //insert the picked children one into each of the list
	GinsertLstNd(ptrChildTemp->ptrChildLst, lstNdOne);
	GsetRect(ptrChildTemp, lstNdOne->tnInfo);
	GinsertLstNd(ptrChildSib->ptrChildLst, lstNdTwo);
	GsetRect(ptrChildSib, lstNdTwo->tnInfo);

	Region rgnNewRectOne = GinitRgnRect(NULL, NULL);
	Region rgnNewRectTwo = GinitRgnRect(NULL, NULL);

	Boolean bIsOne = FALSE;
	Boolean bIsNdOneInComp = FALSE;
	Boolean bIsNdTwoInComp = FALSE;

	int iCnt = 0;

	lstNdTemp = GdeleteLstFirst(ptrChild->ptrChildLst);

    //pick one element from the list of children of the node to be split
	while(lstNdTemp != NULL)
    {   //if one of the nodes has so few entires that all the remaining children are to be assigned set that node to be incomplete
		if(ptrChildTemp->ptrChildLst->uiCnt + ptrChild->ptrChildLst->uiCnt == gminEntries - 1)
			bIsNdOneInComp = TRUE;

		if(ptrChildSib->ptrChildLst->uiCnt + ptrChild->ptrChildLst->uiCnt == gminEntries - 1)
			bIsNdTwoInComp = TRUE;
        //if both nodes are not potentiall incomplete i.e. all the remaining children need not be assigned to it for the node not to underflow
		if(!bIsNdOneInComp && !bIsNdTwoInComp)
        {   dExpOne = -1;
		    dExpTwo = -1;
            //find the area by which the rectangle in each node should be expanded to accomodat the given rectangle
		    GexpansionArea(ptrChildTemp->tnInfo->tdInfo->rgnRect, lstNdTemp->tnInfo, &dExpOne, rgnNewRectOne);
		    GexpansionArea(ptrChildSib->tnInfo->tdInfo->rgnRect, lstNdTemp->tnInfo, &dExpTwo, rgnNewRectTwo);

            //ark the node that requires least expansion to be the list into which the child can be added
		    if(dExpOne < dExpTwo)
			     bIsOne = TRUE;
	        if(dExpOne > dExpTwo)
		         bIsOne = FALSE;
	        if(dExpOne == dExpTwo)
            {    //incase both reequired to be enlarged by same amount pick the one with lower area currently
		         double dAreaOne = Garea(ptrChildTemp->tnInfo->tdInfo->rgnRect);
			     double dAreaTwo = Garea(ptrChildSib->tnInfo->tdInfo->rgnRect);
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
			GinsertLstNd(ptrChildTemp->ptrChildLst, lstNdTemp);
            if(bIsNdOneInComp)
            {   dExpOne = -1;
				GexpansionArea(ptrChildTemp->tnInfo->tdInfo->rgnRect, lstNdTemp->tnInfo, &dExpOne, rgnNewRectOne);
			}
            //expand the rectangle to accomodate new child
			for(iCnt = 0; iCnt < DIMENSION; iCnt++)
            {	ptrChildTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] = rgnNewRectOne->iBottomLeft[iCnt];
				ptrChildTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt] = rgnNewRectOne->iTopRight[iCnt];
			}
		}
		else
        {   //insert the new child
			GinsertLstNd(ptrChildSib->ptrChildLst, lstNdTemp);
            if(bIsNdTwoInComp)
            {   dExpTwo = -1;
				GexpansionArea(ptrChildSib->tnInfo->tdInfo->rgnRect, lstNdTemp->tnInfo, &dExpTwo, rgnNewRectTwo);
			}
            //expand the rectangle to accomodate the new child
			for(iCnt = 0; iCnt < DIMENSION; iCnt++)
            {	ptrChildSib->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] = rgnNewRectTwo->iBottomLeft[iCnt];
				ptrChildSib->tnInfo->tdInfo->rgnRect->iTopRight[iCnt] = rgnNewRectTwo->iTopRight[iCnt];
			}
		}
		lstNdTemp = GdeleteLstFirst(ptrChild->ptrChildLst);
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

void GpickSeeds(GHdrNd ptrChildLst, GLstNd *lstNdChildOne, GLstNd *lstNdChildTwo)
{    
	if(ptrChildLst == NULL)
		return;
	GTreeNode *tnInfoMin = (GTreeNode *)malloc(DIMENSION * sizeof(GTreeNode));
	GTreeNode *tnInfoMax = (GTreeNode *)malloc(DIMENSION * sizeof(GTreeNode));

	GLstNd lstNdTemp = NULL;
    int iCnt = 0;
	Boolean bIsFirst = TRUE;

	Region rgnRectTemp;
	Region rgnRectOut;

	double dNormSep = 0;
	double dMaxNormSep = 0;
	dimension dimMaxNormSep = 0;
	GTreeNode temp;

	switch(ptrChildLst->ptrFirstNd->tnInfo->ndType)
    {	case INTNODE:   lstNdTemp = ptrChildLst->ptrFirstNd;
	                    rgnRectTemp = GinitRgnRect(NULL, NULL);

	                    rgnRectOut = GinitRgnRect(NULL, NULL);
                        while(lstNdTemp != NULL)
                        {     for(iCnt = 0; iCnt < DIMENSION; iCnt++)
                              {     
                              	if(bIsFirst)
                                    {     
				                          rgnRectTemp->iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt];
				                          rgnRectTemp->iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt];
				                          rgnRectOut->iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt];
				                          rgnRectOut->iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt];
				                          tnInfoMin[iCnt] = lstNdTemp->tnInfo;
				                          tnInfoMax[iCnt] = lstNdTemp->tnInfo;
				                          continue;
                                    }
			                        if(lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] > rgnRectTemp->iBottomLeft[iCnt])
                                    {     
				                          rgnRectTemp->iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt];
				                          tnInfoMin[iCnt] = lstNdTemp->tnInfo;
                                    }
			                        if(lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt] < rgnRectTemp->iTopRight[iCnt])
                                    {     
				                          rgnRectTemp->iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight[iCnt];
				                          tnInfoMax[iCnt] = lstNdTemp->tnInfo;
                                    }

			                        if(lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt] < rgnRectOut->iBottomLeft[iCnt])
                                                rgnRectOut->iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft[iCnt];
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

	                     *lstNdChildOne = GdeleteLstElem(ptrChildLst, tnInfoMin[(int)dimMaxNormSep]);

	                     *lstNdChildTwo = GdeleteLstElem(ptrChildLst, tnInfoMax[(int)dimMaxNormSep]);

	                     free(rgnRectTemp->iBottomLeft);
	                     free(rgnRectTemp->iTopRight);
	                     free(rgnRectTemp);
	                     free(rgnRectOut->iBottomLeft);
	                     free(rgnRectOut->iTopRight);
	                     free(rgnRectOut);
                         break;


	case EXTNODE:        
						lstNdTemp = ptrChildLst->ptrFirstNd;
	                    rgnRectTemp = GinitRgnRect(NULL, NULL);
	                    rgnRectOut = GinitRgnRect(NULL, NULL);
                        while(lstNdTemp != NULL)
                        {     for(iCnt = 0; iCnt < DIMENSION; iCnt++)
                              {     if(bIsFirst)
                                    {     
				                          rgnRectTemp->iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->groupElem->bcell->minOriginalBoundary[iCnt];
				                          rgnRectTemp->iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->groupElem->bcell->maxOriginalBoundary[iCnt];
				                          rgnRectOut->iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->groupElem->bcell->minOriginalBoundary[iCnt];
				                          rgnRectOut->iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->groupElem->bcell->maxOriginalBoundary[iCnt];
				                          tnInfoMin[iCnt] = lstNdTemp->tnInfo;
				                          tnInfoMax[iCnt] = lstNdTemp->tnInfo;
				                          continue;
                                    }
			                        if(lstNdTemp->tnInfo->tdInfo->groupElem->bcell->minOriginalBoundary[iCnt] > rgnRectTemp->iBottomLeft[iCnt])
                                    {   
				                          rgnRectTemp->iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->groupElem->bcell->minOriginalBoundary[iCnt];
				                          tnInfoMin[iCnt] = lstNdTemp->tnInfo;
                                    }
			                        if(lstNdTemp->tnInfo->tdInfo->groupElem->bcell->maxOriginalBoundary[iCnt] < rgnRectTemp->iTopRight[iCnt])
                                    {     
				                          rgnRectTemp->iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->groupElem->bcell->maxOriginalBoundary[iCnt];
				                          tnInfoMax[iCnt] = lstNdTemp->tnInfo;
                                    }
			                        if(lstNdTemp->tnInfo->tdInfo->groupElem->bcell->minOriginalBoundary[iCnt] < rgnRectOut->iBottomLeft[iCnt])
                                                rgnRectOut->iBottomLeft[iCnt] = lstNdTemp->tnInfo->tdInfo->groupElem->bcell->minOriginalBoundary[iCnt];
			                        if(lstNdTemp->tnInfo->tdInfo->groupElem->bcell->maxOriginalBoundary[iCnt] > rgnRectOut->iTopRight[iCnt])
                                                rgnRectOut->iTopRight[iCnt] = lstNdTemp->tnInfo->tdInfo->groupElem->bcell->maxOriginalBoundary[iCnt];
                               }
                               lstNdTemp = lstNdTemp->ptrNextNd;
		                       if(bIsFirst)
                                    bIsFirst = FALSE;
                         }	//while
	                     dNormSep = 0;
	                     dMaxNormSep = 0;
	                     dimMaxNormSep = 0;
	                     for(iCnt = 0; iCnt < DIMENSION; iCnt++)
                         {     
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
                         		rgnRectTemp->iTopRight[(int)dimMaxNormSep] = lstNdTemp->tnInfo->tdInfo->groupElem->bcell->maxOriginalBoundary[(int)dimMaxNormSep];
                         		 tnInfoMax[(int)dimMaxNormSep] = lstNdTemp->tnInfo;
                         		lstNdTemp = lstNdTemp->ptrNextNd;
                         	}
                         	else
                         	{
                         		rgnRectTemp->iTopRight[(int)dimMaxNormSep] = lstNdTemp->ptrNextNd->tnInfo->tdInfo->groupElem->bcell->maxOriginalBoundary[(int)dimMaxNormSep];
                         		tnInfoMax[(int)dimMaxNormSep] = lstNdTemp->ptrNextNd->tnInfo;
                         		lstNdTemp = lstNdTemp->ptrNextNd->ptrNextNd;
                         	}
                         	while(lstNdTemp != NULL)
                         	{
                         		if(lstNdTemp->tnInfo->tdInfo->groupElem->bcell->maxOriginalBoundary[(int)dimMaxNormSep] < rgnRectTemp->iTopRight[(int)dimMaxNormSep] && temp != lstNdTemp->tnInfo)
                         		{
                         			rgnRectTemp->iTopRight[(int)dimMaxNormSep] = lstNdTemp->tnInfo->tdInfo->groupElem->bcell->maxOriginalBoundary[(int)dimMaxNormSep];
				                          tnInfoMax[(int)dimMaxNormSep] = lstNdTemp->tnInfo;
                         		}
                         			lstNdTemp = lstNdTemp->ptrNextNd;
                         	} 
                         }

                         // remove the node pointed by tnInfoMin at tnInfoMax at dMaxNormSep and assign them to be the two nodes that are picked for split
	                     *lstNdChildOne = GdeleteLstElem(ptrChildLst, tnInfoMin[(int)dimMaxNormSep]);
	                     *lstNdChildTwo = GdeleteLstElem(ptrChildLst, tnInfoMax[(int)dimMaxNormSep]);

	                     free(rgnRectTemp->iBottomLeft);
	                     free(rgnRectTemp->iTopRight);
	                     free(rgnRectTemp);
	                     free(rgnRectOut->iBottomLeft);
	                     free(rgnRectOut->iTopRight);
	                     free(rgnRectOut);
                         break;
	}//switch

	free(tnInfoMin);
	free(tnInfoMax);

	return;
}

long double distance(DataPoint point1, DataPoint point2)
{
	long double sum = 0;
	int i;
	for(i=0;i<DIMENSION;i++)
		sum += ( (point1[i]- point2[i])*(point1[i] - point2[i]) );
	
	return sqrt(sum);
}

int GisEpsDistant(Region rgnRect, Group groupElem, double eps)
{   
	int iCnt = 0;
	bool assignflag = false;

	int bIsContains = 0;
	double* point = (Dimension)malloc(sizeof(double)*DIMENSION);

	for(iCnt = 0; iCnt < DIMENSION; iCnt++)
    {   
    	point[iCnt] = (rgnRect->iBottomLeft[iCnt] + rgnRect->iTopRight[iCnt])/(double)2;
	}

	if(distance(point, groupElem->master) <= eps)bIsContains = 1;
	else if(distance(point, groupElem->master) < 2*eps)bIsContains = -1;

	free(point);
	return bIsContains;
}

int GcanAssign(GHdrNd hdrNdTree, Region rgnRect, Region x2rgnRect, double eps)
{
	if(GisLstEmpty(hdrNdTree) || rgnRect == NULL || rgnRect->iBottomLeft == NULL || rgnRect->iTopRight == NULL)
		return 0;
 	
 	unsigned int uiRecCnt = 0;
	int iCnt = 0;
    double t;
    static int flag = 0;
	GLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
	int aflag = -2;
	bool eps2 = 0;

	while(lstNdTemp != NULL)
    {	
    	int assignflag;
		switch(lstNdTemp->tnInfo->ndType)
        {   

        	case INTNODE:   //incase of internal node if the node is overlapping with search rectangle descend into the node and add to the count the number of new records found if any
				            if(GisOverLap(lstNdTemp->tnInfo->tdInfo->rgnRect, x2rgnRect))
				            {
				            	assignflag = GcanAssign(lstNdTemp->ptrChildLst, rgnRect, x2rgnRect, eps);
				            	if(assignflag==1) return assignflag;
				            	if(assignflag==-1) eps2 = -1;
				          	}              
                            break;
            case EXTNODE:   
            				aflag = GisEpsDistant(x2rgnRect, lstNdTemp->tnInfo->tdInfo->groupElem, eps);
            				if(aflag == 1)return 1;
            				           				    				
	                        break;

		}//switch
		lstNdTemp = lstNdTemp->ptrNextNd;
	}  

	return eps2;
}
Group GfindCell(GHdrNd hdrNdTree, Region rgnRect)
{   
	//finds the record in the given search rectangle and returns  number of records found
    //if datapoint is passed in iData finds the number of records in eps neighborhood
	if(GisLstEmpty(hdrNdTree) || rgnRect == NULL || rgnRect->iBottomLeft == NULL || rgnRect->iTopRight == NULL)
		return NULL;
 	
 	unsigned int uiRecCnt = 0;
	int iCnt = 0;
    double t;
    static int flag = 0;
	GLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
	Group currGroup = NULL;
	while(lstNdTemp != NULL)
    {	

		switch(lstNdTemp->tnInfo->ndType)
        {   case INTNODE:   //incase of internal node if the node is overlapping with search rectangle descend into the node and add to the count the number of new records found if any
				            if(GisOverLap(lstNdTemp->tnInfo->tdInfo->rgnRect, rgnRect))
				            {
				            	currGroup = GfindCell(lstNdTemp->ptrChildLst, rgnRect);

				            	if(currGroup!=NULL)
				            	{
				            		return currGroup;
				            	}
				            }
                            break;
            case EXTNODE:   
            				if(GisEpsDistant(rgnRect, lstNdTemp->tnInfo->tdInfo->groupElem, EPS)==1)
            				{
            					return lstNdTemp->tnInfo->tdInfo->groupElem;
            				}      			            				
	                        break;

		}//switch
		lstNdTemp = lstNdTemp->ptrNextNd;
	}  
	return NULL;
}

void GprintTree(GHdrNd hdrNdTree)
{   
	//travers along the tree and print the tree
	GLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
	int iCnt = 0;
	static int iIndent = 0;
	iIndent++;

	while(lstNdTemp != NULL)
    {   
    	for(iCnt = 0; iCnt < iIndent; iCnt++)
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
			GprintTree(lstNdTemp->ptrChildLst);
		}
		else
        {   printf("Cell:\n");
    		
    		BCell currBCell = lstNdTemp->tnInfo->tdInfo->groupElem->bcell;
    		
			if(currBCell->cellDataList->count!=0)
			{
				printCell(currBCell);								
			}
			else
			{
				printf("There are no points in this Cell\n");
			}

			printf("End of this Cell:\n ");    				
    		
		}
		lstNdTemp = lstNdTemp->ptrNextNd;
	}
	iIndent--;
	return;
}

void GprintRegion(Region region)
{
    int i=0;
    printf("Printing Region:\n");

    printf("iBottomLeft is: ");
    for (i=0;i<DIMENSION;i++)
    {
        printf("%lf ", region->iBottomLeft[i]);
        
    }
    
    printf("iTopRight is: ");
    for (i=0;i<DIMENSION;i++)
    {
        printf("%lf ", region->iTopRight[i]);
        
    }
    printf("\n");

    return;
    
}

void GgetNeighborHood(GHdrNd hdrNdTree, Data currDataPoint)
{
	// Region pointRgn = getEpsExtendedRegionPoint(currDataPoint,EPS);

	// BCellListHd pointNbhCells = GgetCellsInRegion(hdrNdTree,pointRgn,NULL);
	
	// Data temp = (Data) malloc(sizeof(*temp));
	// *temp = *currDataPoint;
	// int i=0; int val=0;

	// currDataPoint->neighbors = RinitNbHdr();

	// BCellListNode currCellNode = pointNbhCells->first;

	// BCell currBCell;

	// RHdrNd currAuxRTree;

	// while(currCellNode!=NULL)
	// {
	// 	currBCell = currCellNode->bCellElem;
	// 	currAuxRTree = currBCell->auxCellRTree;

	// 	val = RgetNeighborHood(currAuxRTree,temp,0);

	// 	//nbh of currDataPoint is to be appended to the temp
	// 	appendNbh(currDataPoint,temp);
	// 	currCellNode = currCellNode->next;
	// }

	// free(pointRgn->iBottomLeft);
 //    free(pointRgn->iTopRight);
 //    free(pointRgn);

 //    free(temp);

	// freeCellsList(pointNbhCells);

	// return;

}

void appendNbh(Data currDataPoint, Data temp)
{
	if(currDataPoint->neighbors->nbFirst==NULL)
	{
		free(currDataPoint->neighbors);
		currDataPoint->neighbors=temp->neighbors;
		
		return;
	}
	else if (temp->neighbors->nbFirst==NULL)
	{
		free(temp->neighbors);
		return;
	}
	else
	{
		currDataPoint->neighbors->nbLast->nbNext=temp->neighbors->nbFirst;
		currDataPoint->neighbors->nbLast = temp->neighbors->nbLast;
		currDataPoint->neighbors->nbhCnt+=temp->neighbors->nbhCnt;
		free(temp->neighbors);
		return;

	}

	return;

}


BCellListHd GgetCellsInRegion(GHdrNd hdrNdTree, Region regRect, Region cellRgnRect)
{
	// retrieve all the cells overlapping with this region
	BCellListHd cellsList = initBCellListHd();

	int noOfCells = GfindOverlappingCells(hdrNdTree,regRect,cellRgnRect,cellsList);
	return cellsList;
}

// this function excludes the border points located in remote cells falling out of the eps extended region when eps = cellsize
BCellListHd GgetCellsInRegion2(GHdrNd hdrNdTree, Region regRect, Region cellRgnRect)
{
	// retrieve all the cells overlapping with this region
	BCellListHd cellsList = initBCellListHd();

	int noOfCells = GfindOverlappingCells2(hdrNdTree,regRect,cellRgnRect,cellsList);
	return cellsList;
}
BCellListHd GgetCellsInTotalRegion2(GHdrNd hdrNdTree, Region regRect, Region cellRgnRect)
{
	// retrieve all the cells overlapping with this region
	BCellListHd cellsList = initBCellListHd();

	int noOfCells = GfindTotalOverlappingCells2(hdrNdTree,regRect,cellRgnRect,cellsList);

	return cellsList;
}


int GgetCellsInRegion3(GHdrNd hdrNdTree, Region regRect, Region cellRgnRect)
{
	// retrieve all the cells overlapping with this region
	BCellListHd cellsList2 = initBCellListHd();

	int noOfCells = GfindOverlappingCellsTotalOverlap(hdrNdTree,regRect,cellRgnRect,cellsList2);

	int count = 0;
	BCellListNode currCellNode = cellsList2->first;
	Region currCellRgn; Data currData; int i =0; Boolean flag = FALSE; int j=0;
	while(currCellNode!=NULL)
	{
		// check if the cell is within the boundary of immediate region
		currCellRgn = createRegionofCell(currCellNode->bCellElem);

		if(GisOverLap2(regRect,currCellRgn))
		{
			count = count + (currCellNode->bCellElem->cellDataList)->count;

		}
		else
		{
			//check for all the points in this cell if they are on the boundary of the immediate region
			// and then add the count

			CellData currCellData = currCellNode->bCellElem->cellDataList->first;
			
			while(currCellData!=NULL)
			{
				currData = currCellData->data;
				flag = FALSE;
				for(i=0;i<DIMENSION;i++)
				{
					if(currData->iData[i] > regRect->iTopRight[i])
					{
						flag = TRUE;
						break;
					}
				}

				if(flag == FALSE)
				{
					count++;
				}

				currCellData = currCellData->next;
			}
		}
		free(currCellRgn->iBottomLeft);
		free(currCellRgn->iTopRight);
		free(currCellRgn);

		currCellNode = currCellNode->next;
		
	}
	freeCellsList(cellsList2);
	return count;
}

unsigned int GfindOverlappingCells(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList)
{
	// recursively adds the ovelapping cells to the cellsList

	if(GisLstEmpty(hdrNdTree) || rgnRect == NULL || rgnRect->iBottomLeft == NULL || rgnRect->iTopRight == NULL)
		return 0;

	Region currCellRgn;
	unsigned int uiRecCnt = 0;
	int iCnt = 0;
    double t;
    static int flag = 0;
    int k = 0, currflag=0;
	GLstNd lstNdTemp = hdrNdTree->ptrFirstNd;

	while(lstNdTemp != NULL)
    {	

		switch(lstNdTemp->tnInfo->ndType)
        {   case INTNODE:   //incase of internal node if the node is over lapping with search rectangle descend into the node and add to the count the number of new records found if any
				            if(GisOverLap(lstNdTemp->tnInfo->tdInfo->rgnRect, rgnRect)==TRUE)
                                   uiRecCnt += GfindOverlappingCells(lstNdTemp->ptrChildLst, rgnRect, cellRgnRect, cellsList);

                            break;
            case EXTNODE:  
            				currCellRgn = createRegionofCell(lstNdTemp->tnInfo->tdInfo->groupElem->bcell);
            				if(GisOverLap(currCellRgn, rgnRect)==TRUE)
            				{
            						uiRecCnt++;
            						insertBCellIntoBCellList(lstNdTemp->tnInfo->tdInfo->groupElem->bcell,cellsList);
            				}

            				free(currCellRgn->iBottomLeft);
        					free(currCellRgn->iTopRight);
        					free(currCellRgn);
                                                                  
	                        break;

		}//switch
		lstNdTemp = lstNdTemp->ptrNextNd;
	}  
	return uiRecCnt;

}

Boolean GisOverLapTotal_MN(Region rgnRectOne, Region rgnRectTwo, BCell tempCell, int myrank)
{
	// returns TRUE ONLY in case of TOTAL overlap 
	int iDim=0;

	for(iDim=0; iDim < DIMENSION; iDim++)
	{
		if(rgnRectOne->iBottomLeft[iDim] < rgnRectTwo->iBottomLeft[iDim] || rgnRectOne->iTopRight[iDim] > rgnRectTwo->iTopRight[iDim])
		{
			return FALSE;
		}
	}
	return TRUE;
}
unsigned int GfindOverlappingCellsTotalOverlap_MN(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList, BCell tempCell, int myrank)
{
	// recursively adds the ovelapping cells to the cellsList

	if(GisLstEmpty(hdrNdTree) || rgnRect == NULL || rgnRect->iBottomLeft == NULL || rgnRect->iTopRight == NULL)
		return 0;

	Region currCellRgn;
	unsigned int uiRecCnt = 0;
	int iCnt = 0;
    double t;
    static int flag = 0;
    int countPoints=0;
    int k = 0, currflag=0;
	GLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
	// int nodesVisited;
	while(lstNdTemp != NULL)
    {	
    	//nodesVisited++;
		switch(lstNdTemp->tnInfo->ndType)
        {   case INTNODE:   //incase of internal node if the node is over lapping with search rectangle descend into the node and add to the count the number of new records found if any
				            if(GisOverLap(lstNdTemp->tnInfo->tdInfo->rgnRect, rgnRect)==TRUE)
                                   uiRecCnt += GfindOverlappingCellsTotalOverlap_MN(lstNdTemp->ptrChildLst, rgnRect, cellRgnRect, cellsList, tempCell,myrank);

                            break;
            case EXTNODE:  
            				currCellRgn = createRegionofCell(lstNdTemp->tnInfo->tdInfo->groupElem->bcell);
            				if(GisOverLapTotal_MN(currCellRgn, rgnRect, tempCell,myrank)==TRUE)
            				{
            					// Region currCellRgn = createRegionofCell(lstNdTemp->tnInfo->tdInfo->bCellElem);

            						uiRecCnt++;

            						insertBCellIntoBCellList(lstNdTemp->tnInfo->tdInfo->groupElem->bcell,cellsList);
            						countPoints += lstNdTemp->tnInfo->tdInfo->groupElem->bcell->cellDataList->count;   					
            						


            				}
            				else if(GisOverLap(currCellRgn, rgnRect)==TRUE)
            				{
            					int i;
            					Boolean confirm;
            					CellData itr= lstNdTemp->tnInfo->tdInfo->groupElem->bcell->cellDataList->first;
            					while(itr !=NULL)
            					{
            						confirm=TRUE;
            						for(i=0;i<DIMENSION;i++)
            						{
            							if(itr->data->iData[0] < currCellRgn->iBottomLeft[i] || itr->data->iData[0] > currCellRgn->iTopRight[i])
            								confirm=FALSE;
            						}
            						if(confirm == TRUE)
            							countPoints++;

            						itr=itr->next;
            					}



            				}
            				free(currCellRgn->iBottomLeft);
        					free(currCellRgn->iTopRight);
        					free(currCellRgn);
                                                          
	                        break;

		}//switch
		lstNdTemp = lstNdTemp->ptrNextNd;
	}  
	return countPoints;

}
int GgetCellsInRegion_MN(GHdrNd hdrNdTree, Region regRect, Region cellRgnRect, BCell tempCell, int myrank)
{
	// retrieve all the cells overlapping with this region
	BCellListHd cellsList2 = initBCellListHd();

	int noOfCells = GfindOverlappingCellsTotalOverlap_MN(hdrNdTree,regRect,cellRgnRect,cellsList2, tempCell, myrank);

	int count = 0;
	BCellListNode currCellNode = cellsList2->first;
	Region currCellRgn; Data currData; int i =0; Boolean flag = FALSE; int j=0;
	while(currCellNode!=NULL)
	{
		// check if the cell is within the boundary of immediate region
		currCellRgn = createRegionofCell(currCellNode->bCellElem);
		if(GisOverLap2(regRect,currCellRgn))
		{
			count = count + (currCellNode->bCellElem->cellDataList)->count;

		}
		else
		{
			//check for all the points in this cell if they are on the boundary of the immediate region
			// and then add the count

			CellData currCellData = currCellNode->bCellElem->cellDataList->first;
			
			while(currCellData!=NULL)
			{
				currData = currCellData->data;
				flag = FALSE;
				for(i=0;i<DIMENSION;i++)
				{
					if(currData->iData[i] > regRect->iTopRight[i])
					{
						flag = TRUE;
						break;
					}
				}

				if(flag == FALSE)
				{
					count++;
				}

				currCellData = currCellData->next;
			}
		}
		free(currCellRgn->iBottomLeft);
		free(currCellRgn->iTopRight);
		free(currCellRgn);

		currCellNode = currCellNode->next;
		
	}
	freeCellsList(cellsList2);
	return count;
}
unsigned int GfindOverlappingCellsTotalOverlap(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList)
{
	// recursively adds the ovelapping cells to the cellsList

	if(GisLstEmpty(hdrNdTree) || rgnRect == NULL || rgnRect->iBottomLeft == NULL || rgnRect->iTopRight == NULL)
		return 0;

	Region currCellRgn;
	unsigned int uiRecCnt = 0;
	int iCnt = 0;
    double t;
    static int flag = 0;
    int countPoints=0;
    int k = 0, currflag=0;
	GLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
	// int nodesVisited;
	while(lstNdTemp != NULL)
    {	
    	//nodesVisited++;
		switch(lstNdTemp->tnInfo->ndType)
        {   case INTNODE:   //incase of internal node if the node is over lapping with search rectangle descend into the node and add to the count the number of new records found if any
				            if(GisOverLap(lstNdTemp->tnInfo->tdInfo->rgnRect, rgnRect)==TRUE)
                                   uiRecCnt += GfindOverlappingCellsTotalOverlap(lstNdTemp->ptrChildLst, rgnRect, cellRgnRect, cellsList);
                            break;
            case EXTNODE:  
            				currCellRgn = createRegionofCell(lstNdTemp->tnInfo->tdInfo->groupElem->bcell);
            				if(GisOverLapTotal(currCellRgn, rgnRect)==TRUE)
            				{
 
            						uiRecCnt++;
            						insertBCellIntoBCellList(lstNdTemp->tnInfo->tdInfo->groupElem->bcell,cellsList);
            						countPoints += lstNdTemp->tnInfo->tdInfo->groupElem->bcell->cellDataList->count;
            				}
            				else if(GisOverLap(currCellRgn, rgnRect)==TRUE)
            				{
            					int i;
            					Boolean confirm;
            					CellData itr= lstNdTemp->tnInfo->tdInfo->groupElem->bcell->cellDataList->first;
            					while(itr !=NULL)
            					{
            						confirm=TRUE;
            						for(i=0;i<DIMENSION;i++)
            						{
            							if(itr->data->iData[0] < currCellRgn->iBottomLeft[i] || itr->data->iData[0] > currCellRgn->iTopRight[i])
            								confirm=FALSE;
            						}
            						if(confirm == TRUE)
            							countPoints++;

            						itr=itr->next;
            					}



            				}
            				free(currCellRgn->iBottomLeft);
        					free(currCellRgn->iTopRight);
        					free(currCellRgn);
                                                                
	                        break;

		}//switch
		lstNdTemp = lstNdTemp->ptrNextNd;
	}  
	return countPoints;

}

unsigned int GfindOverlappingCellsOptimized(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList, Region halfEpsExtended)//dont give partial boundary cells
{
	// recursively adds the ovelapping cells to the cellsList

	if(GisLstEmpty(hdrNdTree) || rgnRect == NULL || rgnRect->iBottomLeft == NULL || rgnRect->iTopRight == NULL)
		return 0;

	Region currCellRgn;
	unsigned int uiRecCnt = 0;
	int iCnt = 0;
    double t;
    static int flag = 0;
    int k = 0, currflag=0;
	GLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
	// int nodesVisited;
	while(lstNdTemp != NULL)
    {	
    	//nodesVisited++;
		switch(lstNdTemp->tnInfo->ndType)
        {   case INTNODE:   //incase of internal node if the node is over lapping with search rectangle descend into the node and add to the count the number of new records found if any
				            if(GisOverLap2(lstNdTemp->tnInfo->tdInfo->rgnRect, rgnRect)==TRUE)
                                   uiRecCnt += GfindOverlappingCellsOptimized(lstNdTemp->ptrChildLst, rgnRect, cellRgnRect, cellsList,halfEpsExtended);
                            break;
            case EXTNODE:  
            				currCellRgn = createRegionofCell(lstNdTemp->tnInfo->tdInfo->groupElem->bcell);
            				if(GisOverLap2(currCellRgn, halfEpsExtended)==TRUE)
            				{
        						uiRecCnt++;
        						insertBCellIntoBCellList(lstNdTemp->tnInfo->tdInfo->groupElem->bcell,cellsList);          					
            				}
            				else if(GisOverLap2(currCellRgn, rgnRect)==TRUE && lstNdTemp->tnInfo->tdInfo->groupElem->bcell->cellDataList->cellType==SPARSE)
            				{
        						uiRecCnt++;
        						insertBCellIntoBCellList(lstNdTemp->tnInfo->tdInfo->groupElem->bcell,cellsList);          					
            				}

            				free(currCellRgn->iBottomLeft);
        					free(currCellRgn->iTopRight);
        					free(currCellRgn);
                                                                  
	                        break;

		}//switch
		lstNdTemp = lstNdTemp->ptrNextNd;
	}  
	return uiRecCnt;

}
BCellListHd GgetCellsInRegionOptimized(GHdrNd hdrNdTree, Region regRect, Region cellRgnRect, Region halfEpsExtended)
{
	// retrieve all the cells overlapping with this region
	BCellListHd cellsList = initBCellListHd();
	int noOfCells = GfindOverlappingCellsOptimized(hdrNdTree,regRect,cellRgnRect,cellsList,halfEpsExtended);
	return cellsList;
}

unsigned int GfindOverlappingCells2(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList)//dont give partial boundary cells
{
	// recursively adds the ovelapping cells to the cellsList

	if(GisLstEmpty(hdrNdTree) || rgnRect == NULL || rgnRect->iBottomLeft == NULL || rgnRect->iTopRight == NULL)
		return 0;

	Region currCellRgn;
	unsigned int uiRecCnt = 0;
	int iCnt = 0;
    double t;
    static int flag = 0;
    int k = 0, currflag=0;
	GLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
	// int nodesVisited;
	while(lstNdTemp != NULL)
    {	
    	//nodesVisited++;
		switch(lstNdTemp->tnInfo->ndType)
        {   case INTNODE:   //incase of internal node if the node is over lapping with search rectangle descend into the node and add to the count the number of new records found if any
				            if(GisOverLap2(lstNdTemp->tnInfo->tdInfo->rgnRect, rgnRect)==TRUE)
                                   uiRecCnt += GfindOverlappingCells2(lstNdTemp->ptrChildLst, rgnRect, cellRgnRect, cellsList);
                            break;
            case EXTNODE:  
            				currCellRgn = createRegionofCell(lstNdTemp->tnInfo->tdInfo->groupElem->bcell);
            				if(GisOverLap2(currCellRgn, rgnRect)==TRUE)
            				{

            						uiRecCnt++;
            						insertBCellIntoBCellList(lstNdTemp->tnInfo->tdInfo->groupElem->bcell,cellsList);
            				}

            				free(currCellRgn->iBottomLeft);
        					free(currCellRgn->iTopRight);
        					free(currCellRgn);
                                                                  
	                        break;

		}//switch
		lstNdTemp = lstNdTemp->ptrNextNd;
	}  
	return uiRecCnt;

}
unsigned int GfindTotalOverlappingCells2(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList)//dont give partial boundary cells
{
	// recursively adds the ovelapping cells to the cellsList

	if(GisLstEmpty(hdrNdTree) || rgnRect == NULL || rgnRect->iBottomLeft == NULL || rgnRect->iTopRight == NULL)
		return 0;

	Region currCellRgn;
	unsigned int uiRecCnt = 0;
	int iCnt = 0;
    double t;
    static int flag = 0;
    int k = 0, currflag=0;
	GLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
	// int nodesVisited;
	while(lstNdTemp != NULL)
    {	
    	//nodesVisited++;
		switch(lstNdTemp->tnInfo->ndType)
        {   case INTNODE:   //incase of internal node if the node is over lapping with search rectangle descend into the node and add to the count the number of new records found if any
				            if(GisOverLap2(lstNdTemp->tnInfo->tdInfo->rgnRect, rgnRect)==TRUE)
                                   uiRecCnt += GfindTotalOverlappingCells2(lstNdTemp->ptrChildLst, rgnRect, cellRgnRect, cellsList);
                            break;
            case EXTNODE:  
            				currCellRgn = createRegionofCell(lstNdTemp->tnInfo->tdInfo->groupElem->bcell);
            				//if(GisOverLap2(currCellRgn, rgnRect)==TRUE)
            				if(GisOverLapTotal(currCellRgn, rgnRect)==TRUE)
            				{

            						uiRecCnt++;
            						insertBCellIntoBCellList(lstNdTemp->tnInfo->tdInfo->groupElem->bcell,cellsList);
            				}

            				free(currCellRgn->iBottomLeft);
        					free(currCellRgn->iTopRight);
        					free(currCellRgn);
                                                                  
	                        break;

		}//switch
		lstNdTemp = lstNdTemp->ptrNextNd;
	}  
	return uiRecCnt;

}
unsigned int GfindOverlappingCells3(GHdrNd hdrNdTree, Region rgnRect, Region cellRgnRect, BCellListHd cellsList,BCellListHd cellsList2)
{
	// recursively adds the ovelapping cells to the cellsList

	if(GisLstEmpty(hdrNdTree) || rgnRect == NULL || rgnRect->iBottomLeft == NULL || rgnRect->iTopRight == NULL || cellRgnRect == NULL || cellRgnRect->iBottomLeft == NULL || cellRgnRect->iTopRight == NULL)
		return 0;

	Region currCellRgn;
	unsigned int uiRecCnt = 0;
	int iCnt = 0;
    double t;
    static int flag = 0;
    int k = 0, currflag=0;
	GLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
	// int nodesVisited;
	while(lstNdTemp != NULL)
    {	
    	//nodesVisited++;
		switch(lstNdTemp->tnInfo->ndType)
        {   case INTNODE:   //incase of internal node if the node is over lapping with search rectangle descend into the node and add to the count the number of new records found if any
				            if(GisOverLap2(lstNdTemp->tnInfo->tdInfo->rgnRect, rgnRect)==TRUE)
                                   uiRecCnt += GfindOverlappingCells3(lstNdTemp->ptrChildLst, rgnRect, cellRgnRect, cellsList,cellsList2);
                            break;
            case EXTNODE:  
            				currCellRgn = createRegionofCell(lstNdTemp->tnInfo->tdInfo->groupElem->bcell);
            				if(GisOverLap2(currCellRgn, rgnRect)==TRUE)
            				{
            					uiRecCnt++;
            					insertBCellIntoBCellList(lstNdTemp->tnInfo->tdInfo->groupElem->bcell,cellsList);

            				}

            				if(GisOverLap2(currCellRgn, cellRgnRect)==TRUE)
            				{
            						uiRecCnt++;
            						insertBCellIntoBCellList(lstNdTemp->tnInfo->tdInfo->groupElem->bcell,cellsList2);
            				}

            				free(currCellRgn->iBottomLeft);
        					free(currCellRgn->iTopRight);
        					free(currCellRgn);                               
	                        break;

		}//switch
		lstNdTemp = lstNdTemp->ptrNextNd;
	}  
	return uiRecCnt;

}

GHdrNd populateAuxGridRTree(GroupListHd cellsList, int gMinEntries, int gMaxEntries)
{
	int i = 0;
	GroupListNode currGroupListNode = cellsList->first;
	BCell currBCell;
 	Group currGroup;
	GHdrNd hdrNdTree = GinitHdrNd();
	hdrNdTree->ptrFirstNd = GinitLstNd(GinitIntNd(NULL, NULL));
	hdrNdTree->uiCnt++;
	hdrNdTree->ptrFirstNd->ptrChildLst = GinitHdrNd();

	for(i=0;i<cellsList->count;i++)
	{
		currGroup = currGroupListNode->groupElem;
		hdrNdTree = insertGroupIntoRTree(hdrNdTree,currGroup,gMinEntries,gMaxEntries);
		currGroupListNode = currGroupListNode->next;
	}

	return hdrNdTree;

}

vector<int> findNeighbours(DataHdr dataHdrLst, int point_id, double eps)
{

	vector<int> neighbours;
	Data point = dataHdrLst->dataClstElem + point_id;
	int group_id = point->cellId;
	Group g = groupList[group_id];
	int pointcount = g->bcell->cellDataList->count;
	
	CellData data = g->bcell->cellDataList->first;

	while(data != NULL)
	{
		if((point->id != data->data->id) && distance(point->iData, data->data->iData) <= eps)
		{
			neighbours.push_back(data->data->id);
		}
		data = data->next;
	}

	for(int i=0; i < g->reachable_groups.size(); i++)
	{
		Group currGroup = groupList[g->reachable_groups[i]];
		int pointcount = currGroup->bcell->cellDataList->count;
		
		CellData data = currGroup->bcell->cellDataList->first;

		while(data!=NULL)
		{
			if((point->id != data->data->id) && distance(point->iData, data->data->iData) <= eps)
			{
				neighbours.push_back(data->data->id);
			}
			data = data->next;
		}
	}
	return neighbours;
}

void findReachableGroupsofGroupG(GHdrNd GRTree, Group G)
{
		if(GisLstEmpty(GRTree))
		{
			return ;
		}
 	
	 	unsigned int uiRecCnt = 0;
		int iCnt = 0;
	    double t;
	    static int flag = 0;

		GLstNd lstNdTemp = GRTree->ptrFirstNd;
		Group currGroup = NULL;
		Region x3rgnRect = createCellRegOfDataPoint(G->master, 3*EPS);

		while(lstNdTemp != NULL)
	    {
			switch(lstNdTemp->tnInfo->ndType)
	        {   case INTNODE:   //incase of internal node if the node is overlapping with search rectangle descend into the node and add to the count the number of new records found if any
					            if(GisOverLap(lstNdTemp->tnInfo->tdInfo->rgnRect, x3rgnRect))
					            {
					            	findReachableGroupsofGroupG(lstNdTemp->ptrChildLst, G);
					            }
	                            break;
	            case EXTNODE:   
	            				if(distance(G->master, lstNdTemp->tnInfo->tdInfo->groupElem->master) <= 3*EPS)
	            				{
	            					if(G->id != lstNdTemp->tnInfo->tdInfo->groupElem->id)
	            						G->reachable_groups.push_back(lstNdTemp->tnInfo->tdInfo->groupElem->id);
	            				}      			            				
		                        break;

			}//switch
			lstNdTemp = lstNdTemp->ptrNextNd;
		}

		free(x3rgnRect->iBottomLeft);
		free(x3rgnRect->iTopRight);
		free(x3rgnRect);
}


void findReachableGroupsRTree(GHdrNd GRTree, vector<Group>& groupList)
{
	for(int i=0;i<groupList.size();i++)
	{
		Group g = groupList[i];
		findReachableGroupsofGroupG(GRTree, g);
	}
}

GHdrNd populateGridRTree(DataHdr dataHdrLst, GroupListHd groupList, int gMinEntries, int gMaxEntries)
{
	// for loop to read every data point, and every data point
	// find the coords of the BCell
	// check whether the Bcell exists in the R Tree
	// if yes, then insert this data point in that BCell
	// else create a new B Cell, insert this datapoint in it
	// and insert that BCell into R Tree.

	// first create an GR Tree
	GHdrNd hdrNdTree = GinitHdrNd();
	hdrNdTree->ptrFirstNd = GinitLstNd(GinitIntNd(NULL, NULL));
	hdrNdTree->uiCnt++;
	hdrNdTree->ptrFirstNd->ptrChildLst = GinitHdrNd();
	// for loop to process every point
	int cnt = 0, i;
    //retirve element one by one and insert them into tree invoking create root incase of root split
	vector<bool> unassigned(dataHdrLst->uiCnt, false);
							// bool unassigned[dataHdrLst->uiCnt];

	for(i=0;i<dataHdrLst->uiCnt;i++)
	{
		// construct the region with the help of the point
		unassigned[i] = false;
		Data currDataPoint;
		currDataPoint = dataHdrLst->dataClstElem + i;

		Region currRegion = createCellRegOfPoint(currDataPoint, EPS);
		Region x2currRegion = createCellRegOfPoint(currDataPoint, 2*EPS);

		int canBeAssigned = GcanAssign(hdrNdTree, currRegion, x2currRegion, EPS);

		if(canBeAssigned==1)
		{
			Group currGroup = GfindCell(hdrNdTree,currRegion);

			//Modify the implementation of GfindCell to return BCell only when distance between data point and center is less than or equal to Epsilon
			if(currGroup != NULL){
				
				insertPointintoGroup(currGroup, currDataPoint);
				currDataPoint->cellId=currGroup->id;			
			}
			else{
				printf("\n\n\nGcanAssign is giving WRONG OUTPUT!!!\n\n\n");
			}
		}
		else if(canBeAssigned==0)
		{
			Group newGroup = initGroup(currDataPoint);
			newGroup->bcell = initBCell(currRegion);
			BCell newBCell = newGroup->bcell;
			insertPointintoGroup(newGroup, currDataPoint);
			
			currDataPoint->cellId = newGroup->id;

			hdrNdTree = insertGroupIntoRTree(hdrNdTree, newGroup, gMinEntries, gMaxEntries);
			insertGroupIntoGroupList(newGroup, groupList);

		}			
		else
		{
			unassigned[i] = true;
		}

		free(currRegion->iBottomLeft);
   		free(currRegion->iTopRight);
    	free(currRegion);
    	free(x2currRegion->iBottomLeft);
   		free(x2currRegion->iTopRight);
    	free(x2currRegion);
	}

	for(i=0;i<dataHdrLst->uiCnt;i++)
	{
		Region currRegion = NULL;
		Region x2currRegion = NULL;
		if(unassigned[i])
		{
			Data currDataPoint;
			currDataPoint = dataHdrLst->dataClstElem + i;

			currRegion = createCellRegOfPoint(currDataPoint, EPS);
			x2currRegion = createCellRegOfPoint(currDataPoint, 2*EPS);

			int canBeAssigned = GcanAssign(hdrNdTree, currRegion, x2currRegion, EPS);
			
			if(canBeAssigned==1)
			{
				unassigned[i] = false;
				Group currGroup = GfindCell(hdrNdTree,currRegion);

				//Modify the implementation of GfindCell to return BCell only when distance between data point and center is less than or equal to Epsilon
				if(currGroup != NULL){
					
					insertPointintoGroup(currGroup, currDataPoint);
					currDataPoint->cellId=currGroup->id;			
				}
				else{
					printf("\n\n\n GcanAssign is giving WRONG OUTPUT\n\n\n");
				}
			}
			else
			{
				// createa a new Group and insert it into the R Tree
				// create a BCell and insert it into the R Tree
				// printf("Creating a new BCell and inserting the new point into it\n");
					
				Group newGroup = initGroup(currDataPoint);//(Group)malloc(sizeof(struct group));
				newGroup->bcell = initBCell(currRegion);
				BCell newBCell = newGroup->bcell;

				insertPointintoGroup(newGroup, currDataPoint);
					
				currDataPoint->cellId = newGroup->id;

				hdrNdTree = insertGroupIntoRTree(hdrNdTree, newGroup, gMinEntries, gMaxEntries);
				insertGroupIntoGroupList(newGroup, groupList);
			}			
		}
		if(currRegion != NULL && x2currRegion != NULL){
			free(currRegion->iBottomLeft);
	   		free(currRegion->iTopRight);
	    	free(currRegion);
	    	free(x2currRegion->iBottomLeft);
	   		free(x2currRegion->iTopRight);
	    	free(x2currRegion);	
    	}
	}
	
	return hdrNdTree;
}

void insertPointintoGroup(Group currGroup, Data currDataPoint)
{
	long double dist = distance(currDataPoint->iData, currGroup->master);
	
	if(currGroup->threshold < dist) currGroup->threshold = dist;
	currDataPoint->group_id = currGroup->id;

	insertPointintoCellDataList(currGroup->bcell->cellDataList, currDataPoint);

	BCell currBCell = currGroup->bcell;
	int x,i;

	if(currBCell->cellDataList->count == 1)
	{
		for(x = 0; x < DIMENSION; x++)
		{
			currBCell->minActualBoundary[x] = currDataPoint->iData[x];
			currBCell->maxActualBoundary[x] = currDataPoint->iData[x];

		}		
	}
	else
	{
		for(i = 0; i< DIMENSION; i++)
		{
			if(currBCell->minActualBoundary[i] > currDataPoint->iData[i])
			{
				currBCell->minActualBoundary[i] = currDataPoint->iData[i];
			}
			if(currBCell->maxActualBoundary[i] < currDataPoint->iData[i])
			{
				currBCell->maxActualBoundary[i] = currDataPoint->iData[i];
			}
		}

	}

	return;
	
}

void insertGroupIntoGroupList(Group newGroup, GroupListHd cellsList)
{	
	groupList.push_back(newGroup);

	GroupListNode newGroupListNode = (GroupListNode) malloc(sizeof(*newGroupListNode)); // MODIFIED	

	newGroupListNode->groupElem = newGroup;	

	newGroupListNode->next = cellsList->first;
	cellsList->first = newGroupListNode;
	cellsList->count++;
	return;
}

void insertBCellIntoBCellList(BCell newBCell,BCellListHd cellsList)
{	
	BCellListNode newBCellListNode = (BCellListNode) malloc(sizeof(*newBCellListNode)); // MODIFIED	
	newBCellListNode->bCellElem = newBCell;	

	newBCellListNode->next = cellsList->first;
	cellsList->first = newBCellListNode;
	cellsList->count++;

	return;

}

void insertPointintoCellDataList(CellDataHd currCellDataList, Data currDataPoint)
{
	CellData currCellData = initCellData(currDataPoint);
	
	if(currCellDataList->first==NULL)
	{
		currCellDataList->first = currCellData;
		currCellDataList->count++;

	}
	else
	{
		currCellData->next = currCellDataList->first;
		currCellDataList->first = currCellData;
		currCellDataList->count++;

	}

	return;

}

Region createCellRegOfDataPoint(DataPoint currDataPoint, double eps)
{
	Region tempRegion = (Region) malloc(sizeof(struct region));
	assert(tempRegion!=NULL);

	tempRegion->iBottomLeft = (Dimension)malloc(DIMENSION * sizeof(double));
	assert(tempRegion->iBottomLeft!=NULL);

	tempRegion->iTopRight = (Dimension)malloc(DIMENSION * sizeof(double));
	assert(tempRegion->iTopRight!=NULL);

	int i = 0;

	for(i=0;i<DIMENSION;i++)
	{
		tempRegion->iBottomLeft[i] = currDataPoint[i] - eps;		// * CELLSIZE;
		tempRegion->iTopRight[i] = currDataPoint[i] + eps; 		// * CELLSIZE;
	}

	return tempRegion;
}

Region createCellRegOfPoint(Data currDataPoint, double eps)
{
	Region tempRegion = (Region) malloc(sizeof(struct region));
	assert(tempRegion!=NULL);

	tempRegion->iBottomLeft = (Dimension)malloc(DIMENSION * sizeof(double));
	assert(tempRegion->iBottomLeft!=NULL);

	tempRegion->iTopRight = (Dimension)malloc(DIMENSION * sizeof(double));
	assert(tempRegion->iTopRight!=NULL);

	int i = 0;
	for(i=0;i<DIMENSION;i++)
	{
		tempRegion->iBottomLeft[i] = currDataPoint->iData[i] - eps;		// * CELLSIZE;
		tempRegion->iTopRight[i] = currDataPoint->iData[i] + eps; 		// * CELLSIZE;
	}

	return tempRegion;
}

Region createRegionofCell(BCell bCellElem)
{
	Region tempRegion = (Region) malloc(sizeof(struct region));
	assert(tempRegion!=NULL);

	tempRegion->iBottomLeft = (Dimension)malloc(DIMENSION * sizeof(double));
	assert(tempRegion->iBottomLeft!=NULL);

	tempRegion->iTopRight = (Dimension)malloc(DIMENSION * sizeof(double));
	assert(tempRegion->iTopRight!=NULL);

	int iCnt;
	for(iCnt = 0; iCnt < DIMENSION; iCnt++)
    {   
    	tempRegion->iBottomLeft[iCnt] = bCellElem->minOriginalBoundary[iCnt];
		tempRegion->iTopRight[iCnt] = bCellElem->maxOriginalBoundary[iCnt];
	}

	return tempRegion;
}

GHdrNd insertGroupIntoRTree(GHdrNd hdrNdTree, Group groupNode, int gMinEntries, int gMaxEntries)
{
	GinsertTree(hdrNdTree, GinitExtNd(groupNode), gMinEntries, gMaxEntries);
		
	if(hdrNdTree->uiCnt > 1)
		hdrNdTree = GcreateRoot(hdrNdTree);

	return hdrNdTree;
}



void printMinGridSize()
{
	int i = 0;
	printf("MINGRIDSIZE is: ");
	for(i=0;i<DIMENSION;i++)
	{
		printf("%lf ",MINGRIDSIZE[i]);
	}
	printf("\n");

	return;
}

void printMaxGridSize()
{
	int i = 0;
	printf("MAXGRIDSIZE is: ");
	for(i=0;i<DIMENSION;i++)
	{
		printf("%lf ",MAXGRIDSIZE[i]);
	}
	printf("\n");

	return;

}

void printCellsList(BCellListHd cellsList)
{
	int i=0;
	BCellListNode currCellNode = cellsList->first;
	printf("cell count=%d\n",cellsList->count);
	
	if(cellsList->count != 0)
	{
		while(currCellNode!=NULL)
		{
			BCell currBCellElem = currCellNode->bCellElem;
			if(currBCellElem->cellDataList->count!=0)
			{				
				printf("Printing Cell:\n");
				printCell(currBCellElem);
				printf("End of Printing Cell\n");
			}
			else
			{
				printf("The current cell doesn't have any datapoints in it\n");
			}
			
			currCellNode = currCellNode->next;
		}

	}
	else
	{
		printf("The CellList doesn't have any cells populated\n");
	}	
	
	return;

}

void printNoOfCoreCells(BCellListHd cellsList)
{
	int i=0;
	int coreCells = 0;
	BCellListNode currCellNode = cellsList->first;
	printf("cell count=%d\n",cellsList->count);
	if(cellsList->count != 0)
	{
		while(currCellNode!=NULL)
		{
			BCell currBCellElem = currCellNode->bCellElem;
			if((currBCellElem->cellDataList->count) >= MINPOINTS)
			{
				coreCells++;
			}
			currCellNode = currCellNode->next;
		}
	}
	else
	{
		printf("The CellList doesn't have any cells populated\n");
	}
	
	printf("No of core cells = %d\n", coreCells);

	return;
}

void printCell(BCell bCellElem)
{
	int iCnt=0;
	CellData currDataPt = bCellElem->cellDataList->first;

	for(iCnt=0;iCnt<DIMENSION;iCnt++)
	{
		printf("%lf ", bCellElem->minOriginalBoundary[iCnt]);
	}
	for(iCnt=0;iCnt<DIMENSION;iCnt++)
	{
		printf("%lf ", bCellElem->maxOriginalBoundary[iCnt]);		
	}
	printf("\n");
	
	while(currDataPt!=NULL)
	{
		printf("%d ", currDataPt->data->iNum);
		for(iCnt=0;iCnt<DIMENSION;iCnt++)
		{
			printf("%lf ", currDataPt->data->iData[iCnt]);
		}
		printf("\n");
		currDataPt = currDataPt->next;
	}

	return;

}

void printCellDataList(CellDataHd currCellDataList)
{
	CellData currCellData = currCellDataList->first;
	printf("\n");
	while(currCellData!=NULL)
	{
		printf("%d ",currCellData->data->iNum);
		currCellData=currCellData->next;
	}

	return;
	
}

void printCellData(CellData currCellData)
{
	int i = 0;
	printf("Printing Cell Data: ");
	printf("%d ",currCellData->data->iNum);

	for(i=0;i<DIMENSION;i++)
	{
	}
	printf("\n");
	return;

}

Region getEpsExtendedRegion(Region cellRgnRect, double eps)
{
	Region tempRegion = (Region) malloc(sizeof(struct region));
	assert(tempRegion!=NULL);

	tempRegion->iBottomLeft = (Dimension)malloc(DIMENSION * sizeof(double));
	assert(tempRegion->iBottomLeft!=NULL);

	tempRegion->iTopRight = (Dimension)malloc(DIMENSION * sizeof(double));
	assert(tempRegion->iTopRight!=NULL);

	int iCnt;
	for(iCnt = 0; iCnt < DIMENSION; iCnt++)
    {   
    	tempRegion->iBottomLeft[iCnt] = cellRgnRect->iBottomLeft[iCnt] - eps;
		tempRegion->iTopRight[iCnt] = cellRgnRect->iTopRight[iCnt] + eps;
	}

	return tempRegion;

}
Region getEpsOptimalExtendedRegion(Region cellRgnRect, double eps)
{
	Region tempRegion = (Region) malloc(sizeof(struct region));
	assert(tempRegion!=NULL);

	tempRegion->iBottomLeft = (Dimension)malloc(DIMENSION * sizeof(double));
	assert(tempRegion->iBottomLeft!=NULL);

	tempRegion->iTopRight = (Dimension)malloc(DIMENSION * sizeof(double));
	assert(tempRegion->iTopRight!=NULL);

	int iCnt;
	for(iCnt = 0; iCnt < DIMENSION; iCnt++)
    {   
    	if(iCnt==0)
    		tempRegion->iBottomLeft[iCnt] = cellRgnRect->iBottomLeft[iCnt];
    	else
    		tempRegion->iBottomLeft[iCnt] = cellRgnRect->iBottomLeft[iCnt] - eps;
    	tempRegion->iTopRight[iCnt] = cellRgnRect->iTopRight[iCnt] + eps;
	}

	return tempRegion;

}

Region getEpsExtendedRegionPoint(Data dataPoint, double eps)
{
	Region tempRegion = (Region) malloc(sizeof(struct region));
	assert(tempRegion!=NULL);

	tempRegion->iBottomLeft = (Dimension)malloc(DIMENSION * sizeof(double));
	assert(tempRegion->iBottomLeft!=NULL);

	tempRegion->iTopRight = (Dimension)malloc(DIMENSION * sizeof(double));
	assert(tempRegion->iTopRight!=NULL);

	int iCnt;
	for(iCnt = 0; iCnt < DIMENSION; iCnt++)
    {   
    	tempRegion->iBottomLeft[iCnt] = dataPoint->iData[iCnt] - eps;
		tempRegion->iTopRight[iCnt] = dataPoint->iData[iCnt] + eps;
	}

	return tempRegion;

}

void freeGRTree(GHdrNd hdrNdTree)
{

	if(hdrNdTree == NULL) 
		return;

	else if(hdrNdTree->uiCnt == 0)
	{
		free(hdrNdTree);
		return;
	}
	
	GLstNd lstNdTemp = hdrNdTree->ptrFirstNd;
	GLstNd lstNdNextTemp;
	
	if(lstNdTemp!=NULL)
	{
		while(lstNdTemp != NULL)
    	{
			switch(lstNdTemp->tnInfo->ndType)
        	{   case INTNODE:   
        						lstNdNextTemp = lstNdTemp->ptrNextNd;
				    			free(lstNdTemp->tnInfo->tdInfo->rgnRect->iBottomLeft);
				    			free(lstNdTemp->tnInfo->tdInfo->rgnRect->iTopRight);
				    			free(lstNdTemp->tnInfo->tdInfo->rgnRect);
				    			free(lstNdTemp->tnInfo->tdInfo);
				    			free(lstNdTemp->tnInfo);
				    			freeGRTree(lstNdTemp->ptrChildLst);                            	
                            	free(lstNdTemp);
                            	break;
            	case EXTNODE:   	
					            lstNdNextTemp = lstNdTemp->ptrNextNd;
					            free(lstNdTemp->tnInfo->tdInfo);
					            free(lstNdTemp->tnInfo);
					            freeGRTree(lstNdTemp->ptrChildLst);
                            	free(lstNdTemp);
		                        break;

			}//switch
			
			lstNdTemp = lstNdNextTemp;
		}
	}	
	
	free(hdrNdTree);	
	return;
}

void isCorrectGRTree(GHdrNd hdrNdTree)
{
    printf("Checking correctness\n");
    GLstNd lstNdTemp = hdrNdTree->ptrFirstNd,temp;
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
            isCorrectGRTree(lstNdTemp->ptrChildLst);
        }
        lstNdTemp = lstNdTemp->ptrNextNd;
    }
    return;
}