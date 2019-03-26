#include "GList.h"

GHdrNd GinitHdrNd()
{   GHdrNd HdrNdLst = (GHdrNd)malloc(sizeof(struct GhdrNd));
	assert(HdrNdLst!=NULL);    

	HdrNdLst->uiCnt = 0;
	HdrNdLst->ptrFirstNd = NULL;
	HdrNdLst->ptrParentNd = NULL;
	
	return HdrNdLst;
}

GLstNd GinitLstNd(GTreeNode tnInfo)
{   GLstNd LstNdElem = (GLstNd)malloc(sizeof(struct GlstNd));
	assert(LstNdElem!=NULL);	

	LstNdElem->tnInfo = tnInfo;
	LstNdElem->ptrChildLst = NULL;
	LstNdElem->ptrNextNd = NULL;
	
	return LstNdElem;	
}

void GinsertLstElem(GHdrNd HdrNdLst, GTreeNode tnInfo)
{	GinsertLstNd(HdrNdLst, GinitLstNd(tnInfo));
	return;
}

void GinsertLstNd(GHdrNd HdrNdLst, GLstNd LstNdElem)
{   if(LstNdElem == NULL || HdrNdLst == NULL)
		return;
    LstNdElem->ptrNextNd = HdrNdLst->ptrFirstNd;
	HdrNdLst->ptrFirstNd = LstNdElem;
	HdrNdLst->uiCnt++;
	
	return;
}

Boolean GisLstEmpty(GHdrNd HdrNdLst)
{   return(HdrNdLst->ptrFirstNd == NULL || HdrNdLst->uiCnt == 0) ? TRUE : FALSE;
}

GLstNd GdeleteLstElem(GHdrNd HdrNdLst, GTreeNode tnInfo)
{ if(GisLstEmpty(HdrNdLst))
    return NULL;
  GTreeNode tnInfoTemp = NULL;
  Boolean bIsFound = FALSE;
  GLstNd LstNdTemp = HdrNdLst->ptrFirstNd;
  if(LstNdTemp->tnInfo == tnInfo)
  { 
    	return GdeleteLstFirst(HdrNdLst);
  }

  while(LstNdTemp->ptrNextNd != NULL)
    {   if(LstNdTemp->ptrNextNd->tnInfo == tnInfo)
        {   bIsFound = TRUE;
      // LstNdTemp = deleteNextNd(HdrNdLst, LstNdTemp);
      // break;     
      return GdeleteNextNd(HdrNdLst, LstNdTemp);
    }   
    LstNdTemp = LstNdTemp->ptrNextNd;
  }
  
    //tnInfoTemp = LstNdTemp->tnInfo;
    //LstNdTemp->tnInfo = NULL;
    //free( LstNdTemp );
    //return (bIsFound)?LstNdTemp:NULL;
  return NULL;
}

GLstNd GdeleteLstFirst(GHdrNd HdrNdLst)
{
   if(GisLstEmpty(HdrNdLst))
		return NULL;
	GTreeNode tnInfoTemp = NULL;
	GLstNd lstNdTemp = HdrNdLst->ptrFirstNd;

	HdrNdLst->ptrFirstNd = lstNdTemp->ptrNextNd;
	HdrNdLst->uiCnt--;
	lstNdTemp->ptrNextNd = NULL;

		
	return lstNdTemp;
}

GLstNd GdeleteNextNd(GHdrNd HdrNdLst, GLstNd LstNdElem)
{
   if(GisLstEmpty(HdrNdLst) || LstNdElem == NULL || LstNdElem->ptrNextNd == NULL)	
		return NULL;

	GLstNd LstNdTemp = LstNdElem->ptrNextNd;
    LstNdElem->ptrNextNd = LstNdTemp->ptrNextNd;
	LstNdTemp->ptrNextNd = NULL;	
	HdrNdLst->uiCnt--;

	return LstNdTemp;
}
