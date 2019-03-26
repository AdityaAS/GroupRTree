#include "RList.h"

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

RHdrNd RinitHdrNd()
{   
	RHdrNd HdrNdLst = (RHdrNd)malloc(sizeof(struct RhdrNd));
	assert(HdrNdLst!=NULL);    

	HdrNdLst->uiCnt = 0;
	HdrNdLst->ptrFirstNd = NULL;
	HdrNdLst->ptrParentNd = NULL;
	
	return HdrNdLst;
}

RLstNd RinitLstNd(RTreeNode tnInfo)
{   
	RLstNd LstNdElem = (RLstNd)malloc(sizeof(struct RlstNd));
	assert(LstNdElem!=NULL);	

	LstNdElem->tnInfo = tnInfo;
	LstNdElem->ptrChildLst = NULL;
	LstNdElem->ptrNextNd = NULL;
	
	return LstNdElem;	
}

void RinsertLstElem(RHdrNd HdrNdLst, RTreeNode tnInfo)
{	
	RinsertLstNd(HdrNdLst, RinitLstNd(tnInfo));
	return;
}

void RinsertLstNd(RHdrNd HdrNdLst, RLstNd LstNdElem)
{   
	if(LstNdElem == NULL || HdrNdLst == NULL)
		return;
    LstNdElem->ptrNextNd = HdrNdLst->ptrFirstNd;
	HdrNdLst->ptrFirstNd = LstNdElem;
	HdrNdLst->uiCnt++;
	
	return;
}

Boolean RisLstEmpty(RHdrNd HdrNdLst)
{   
    if(HdrNdLst!=NULL)
    {
        if(HdrNdLst->ptrFirstNd == NULL)
        return TRUE;

        else if(HdrNdLst->uiCnt == 0)
            return TRUE;

        else
            return FALSE;
    
    }
	
	return TRUE;
    
}

RLstNd RdeleteLstElem(RHdrNd HdrNdLst, RTreeNode tnInfo)
{ if(RisLstEmpty(HdrNdLst))
    return NULL;
  RTreeNode tnInfoTemp = NULL;
  Boolean bIsFound = FALSE;
  RLstNd LstNdTemp = HdrNdLst->ptrFirstNd;
  if(LstNdTemp->tnInfo == tnInfo)
    { return RdeleteLstFirst(HdrNdLst);  }

  while(LstNdTemp->ptrNextNd != NULL)
    {   if(LstNdTemp->ptrNextNd->tnInfo == tnInfo)
        {   bIsFound = TRUE;
        
      return RdeleteNextNd(HdrNdLst, LstNdTemp);
    }   
    LstNdTemp = LstNdTemp->ptrNextNd;
  }
  
  return NULL;
}

RLstNd RdeleteLstFirst(RHdrNd HdrNdLst)
{   
	if(RisLstEmpty(HdrNdLst))
		return NULL;
	RTreeNode tnInfoTemp = NULL;
	RLstNd lstNdTemp = HdrNdLst->ptrFirstNd;

	HdrNdLst->ptrFirstNd = lstNdTemp->ptrNextNd;
	HdrNdLst->uiCnt--;
	lstNdTemp->ptrNextNd = NULL;

		
	return lstNdTemp;
}

RLstNd RdeleteNextNd(RHdrNd HdrNdLst, RLstNd LstNdElem)
{   
	if(RisLstEmpty(HdrNdLst) || LstNdElem == NULL || LstNdElem->ptrNextNd == NULL)	
		return NULL;

	RLstNd LstNdTemp = LstNdElem->ptrNextNd;
    LstNdElem->ptrNextNd = LstNdTemp->ptrNextNd;
	LstNdTemp->ptrNextNd = NULL;	
	HdrNdLst->uiCnt--;

	return LstNdTemp;
}


RNbHdr RinitNbHdr()
{
    RNbHdr nbHdrInfo = (RNbHdr)malloc(sizeof(struct RnbHdr));
    if(nbHdrInfo == NULL)
        return NULL;

    nbHdrInfo->nbhCnt = 0;
    nbHdrInfo->nbFirst =NULL;
    nbHdrInfo->nbLast = NULL;

    return nbHdrInfo;
}

Boolean RisNbLstEmpty(RNbHdr nbHdrInfo)
{    return(nbHdrInfo == NULL || nbHdrInfo->nbhCnt == 0) ? TRUE : FALSE;
}
void RinsertNeighbors(Data dataTemp, Data dataClstElem, double dist)
{
    
   RNbHdr nbHdrInfo = dataTemp->neighbors;

    RNB newNb;

    if(RisNbLstEmpty(nbHdrInfo))
    {   newNb = (RNB)malloc(sizeof(struct Rnb));

        newNb->n = dataClstElem->iNum;
        newNb->Gindex=dataClstElem->Gindex;

        newNb->dist = dist;
        newNb->nbNext = NULL;

        nbHdrInfo->nbFirst = newNb;
        nbHdrInfo->nbLast = newNb;
        
    }
    else
    {   newNb = (RNB)malloc(sizeof(struct Rnb));
        newNb->n = dataClstElem->iNum;
        newNb->Gindex=dataClstElem->Gindex;

        newNb->dist = dist;
        newNb->nbNext = NULL;

        nbHdrInfo->nbLast->nbNext = newNb;
        nbHdrInfo->nbLast = newNb;

    }
    dataTemp->neighbors->nbhCnt++;

    return;
   

}

void RprintNbhLst(RNbHdr nbHdrInfo, DataHdr dataList1)
{   
    if(RisNbLstEmpty(nbHdrInfo))
        return;
    //printf("++++++++++CALL TO PRINTing THE representative  LIST +++++++++++++++++++\n");
    int iCnt = 0;
    RNB nbTemp = nbHdrInfo->nbFirst;

    while(nbTemp != NULL)
    {   //printf("++++++++++CALL TO PRINT THE LIST__________+++++++++++++++++++\n");
        printf("\nN = %d \t dist = %lf\n", nbTemp->n, nbTemp->dist);
        printData(dataList1->dataClstElem + nbTemp->n-1);
        nbTemp = nbTemp->nbNext;
   	}
    return;
}

void RfreeNeighborhood(Data dataTemp)
{
    if(dataTemp->neighbors!=NULL)
    {
        RNB nb_temp = dataTemp->neighbors->nbFirst,nb_next;
        while(nb_temp!=NULL)
        {
            nb_next=nb_temp->nbNext;
            free(nb_temp);
            nb_temp=nb_next;
        }
        free(dataTemp->neighbors);    
    }
	
	return;

}




