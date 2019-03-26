#include "Data.h"
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

void insertDataLstElem(DataHdr dataHdrInfo, DataPoint iData)
{

    Data dataClstElem = dataHdrInfo->dataClstElem + dataHdrInfo->uiCnt;

    dataClstElem->iData = iData;
    dataHdrInfo->uiCnt++;
    dataClstElem->iNum=dataHdrInfo->uiCnt;
    dataClstElem->neighbors = NULL;
    dataClstElem->isProcessed = FALSE;
    dataClstElem->group_id = -1;
    dataClstElem->core_tag = FALSE;
    dataClstElem->neighbors = NULL;
    dataClstElem->ClusterID=UNDEFINED;
    dataClstElem->cellId=UNDEFINED;
    dataClstElem->isProcessed=FALSE;
    dataClstElem->id=DATAID++;
    return;

}

void insertDataLstElem1(DataHdr dataHdrInfo, DataPoint iData, int index, int serial)
{
    Data dataClstElem = dataHdrInfo->dataClstElem + dataHdrInfo->uiCnt;

   dataClstElem->iData = iData;
    dataHdrInfo->uiCnt++;
    dataClstElem->Gindex=serial;
    dataClstElem->iNum=index;

    dataClstElem->neighbors = NULL;
     dataClstElem->isProcessed = FALSE;
    dataClstElem->core_tag = FALSE;
    dataClstElem->neighbors = NULL;
    dataClstElem->ClusterID=UNDEFINED;
    dataClstElem->cellId=UNDEFINED;
    dataClstElem->isProcessed=FALSE;

    return;
}



DataHdr initDataHdr(int size)
{   DataHdr dataHdrInfo = (DataHdr)malloc(sizeof(struct dataHdr));
    assert(dataHdrInfo!=NULL);

    dataHdrInfo->uiCnt = 0;
    dataHdrInfo->dataClstElem=(Data)malloc(sizeof(struct data)*size);
    assert(dataHdrInfo->dataClstElem!=NULL);
    return dataHdrInfo;

}

Data initData(DataPoint iData)
{   
 if(iData == NULL)
        return NULL;

    Data dataClstElem = (Data)malloc(sizeof(struct data));
    assert(dataClstElem!=NULL);
    if(dataClstElem == NULL)
        return NULL;
    dataClstElem->iData = iData;
    dataClstElem->iNum=0;
    dataClstElem->group_id = -1;
    dataClstElem->isProcessed = FALSE;
    dataClstElem->core_tag = FALSE;
    dataClstElem->neighbors = NULL;
    dataClstElem->ClusterID=UNDEFINED;
    dataClstElem->cellId=UNDEFINED;
    dataClstElem->isProcessed=FALSE;
    return dataClstElem;
}

void freeRNbHdr(RNbHdr neighbor)
{
    RNB last = neighbor->nbLast;
    RNB temp = neighbor->nbFirst;
 
    while(temp!=last)
    {
        RNB freed = temp;
        temp = temp->nbNext;

        free(freed);
    }
    free(temp);
    free(neighbor);
    return;
}
void deleteData(DataHdr dataList)
{
    int cnt = dataList->uiCnt;
    for(int i=0;i<cnt;i++)
    {
        Data data = dataList->dataClstElem + i;
        free(data->iData);
    }
    free(dataList->dataClstElem);
    return;
}

DataHdr readData(char *strFileName)
{   
    // read data from the file
    int iCnt = 0;
    int iDim = 0;
    int count=1;
    
    if(strFileName == NULL)
        return NULL;

    FILE *filP = fopen(strFileName, "r");

    if(filP == NULL)
        return NULL;

    DataPoint dataPointTemp = NULL;
    fscanf(filP, "%d", &noOfPoints);
    fscanf(filP, "%d", &DIMENSION);

    MINGRIDSIZE = (double *) malloc(DIMENSION*sizeof(double));
    MAXGRIDSIZE = (double *) malloc(DIMENSION*sizeof(double));
    
    DataHdr dataList = initDataHdr(noOfPoints);

    if(dataList == NULL)
        return NULL;

    dataPointTemp = (DataPoint)malloc(DIMENSION *sizeof(dataPoint));


        for(iDim = 0; iDim < DIMENSION - 1; iDim++)
        {
            fscanf(filP, "%lf"DELIM, &dataPointTemp[iDim]);
            MINGRIDSIZE[iDim] = dataPointTemp[iDim];
            MAXGRIDSIZE[iDim] = dataPointTemp[iDim];
        }            

        fscanf(filP, "%lf\n", &dataPointTemp[iDim]);
        MINGRIDSIZE[iDim] = dataPointTemp[iDim];
        MAXGRIDSIZE[iDim] = dataPointTemp[iDim];
        insertDataLstElem(dataList, dataPointTemp);

        
    while(!feof(filP))
    {   dataPointTemp = (DataPoint)malloc(DIMENSION *sizeof(dataPoint));

        for(iDim = 0; iDim < DIMENSION - 1; iDim++)
        {
            fscanf(filP, "%lf"DELIM, &dataPointTemp[iDim]);
            if(MINGRIDSIZE[iDim] > dataPointTemp[iDim])
                MINGRIDSIZE[iDim] = dataPointTemp[iDim];
            if(MAXGRIDSIZE[iDim] < dataPointTemp[iDim])
                MAXGRIDSIZE[iDim] = dataPointTemp[iDim];
        }
            

        fscanf(filP, "%lf\n", &dataPointTemp[iDim]);
        if(MINGRIDSIZE[iDim] > dataPointTemp[iDim])
                MINGRIDSIZE[iDim] = dataPointTemp[iDim];
            if(MAXGRIDSIZE[iDim] < dataPointTemp[iDim])
                MAXGRIDSIZE[iDim] = dataPointTemp[iDim];
        insertDataLstElem(dataList, dataPointTemp);
    }
    
    fclose(filP);

    return dataList;

}

int readPoints(char *strFileName){

    // read data from the file
    
    
    int noOfPoints;
    
    if(strFileName == NULL)
        return 0;

    FILE *filP = fopen(strFileName, "r");

    if(filP == NULL)
        return 0;
    //DataPoint *ptrData = (DataPoint *) malloc( iNumber * sizeof(DataPoint) );
    
    fscanf(filP,"%d",&noOfPoints);
    

    fclose(filP);

    return noOfPoints;

}

DataHdr readSampledData(char *strFileName)
{   // read data from the file
    
    int iCnt = 0;
    int iDim = 0;
    int count=1;
    int noOfPoints;
    
    if(strFileName == NULL)
        return NULL;

    FILE *filP = fopen(strFileName, "r");

    if(filP == NULL)
        return NULL;
    //DataPoint *ptrData = (DataPoint *) malloc( iNumber * sizeof(DataPoint) );
    DataPoint dataPointTemp = NULL;
    fscanf(filP,"%d",&noOfPoints);
    fscanf(filP,"%d",&DIMENSION);

    // MINGRIDSIZE = (double *) malloc(DIMENSION*sizeof(double));
    // MAXGRIDSIZE = (double *) malloc(DIMENSION*sizeof(double));

    DataHdr dataList = initDataHdr(noOfPoints);

    if(dataList == NULL)
        return NULL;

    dataPointTemp = (DataPoint)malloc(DIMENSION *sizeof(dataPoint));

        for(iDim = 0; iDim < DIMENSION - 1; iDim++)
        {
            fscanf(filP, "%lf"DELIM, &dataPointTemp[iDim]);
            // MINGRIDSIZE[iDim] = dataPointTemp[iDim];
            // MAXGRIDSIZE[iDim] = dataPointTemp[iDim];
        }            

        fscanf(filP, "%lf\n", &dataPointTemp[iDim]);
        // MINGRIDSIZE[iDim] = dataPointTemp[iDim];
        // MAXGRIDSIZE[iDim] = dataPointTemp[iDim];
        insertDataLstElem(dataList, dataPointTemp);



    while(!feof(filP))
    {   dataPointTemp = (DataPoint)malloc(DIMENSION *sizeof(dataPoint));

        for(iDim = 0; iDim < DIMENSION - 1; iDim++)
        {
            fscanf(filP, "%lf"DELIM, &dataPointTemp[iDim]);
            // if(MINGRIDSIZE[iDim] > dataPointTemp[iDim])
            //     MINGRIDSIZE[iDim] = dataPointTemp[iDim];
            // if(MAXGRIDSIZE[iDim] < dataPointTemp[iDim])
            //     MAXGRIDSIZE[iDim] = dataPointTemp[iDim];
        }
            

        fscanf(filP, "%lf\n", &dataPointTemp[iDim]);
        // if(MINGRIDSIZE[iDim] > dataPointTemp[iDim])
        //         MINGRIDSIZE[iDim] = dataPointTemp[iDim];
        //     if(MAXGRIDSIZE[iDim] < dataPointTemp[iDim])
        //         MAXGRIDSIZE[iDim] = dataPointTemp[iDim];
        insertDataLstElem(dataList, dataPointTemp);
    }

    fclose(filP);

    return dataList;
}

void printDataLst(DataHdr dataHdrInfo)
{    //printf("++++++++++CALL TO PRINT THE LIST__________+++++++++++++++++++\n");
    if(isDataLstEmpty(dataHdrInfo))
        return;
    //printf("++++++++++CALL TO PRINTing THE representative  LIST +++++++++++++++++++\n");
    int iCnt = 0,j;
    //DataNd dataTemp = dataHdrInfo->dataFirstNd;
    for(j=0;j<dataHdrInfo->uiCnt;j++)
    {
        printf("RecNum %d\t",dataHdrInfo->dataClstElem[j].iNum);
        
        for(iCnt = 0; iCnt < DIMENSION; iCnt++ )
            printf("%lf ", dataHdrInfo->dataClstElem[j].iData[iCnt] );
            
        //printf("Core Tag %d\t",dataHdrInfo->dataClstElem[j].core_tag);
        //printf("Cluster %d",dataHdrInfo->dataClstElem[j].dataCluster);
        printf("\n");

    }

  return;
}


void printDataLstFile(DataHdr dataHdrInfo,char * strFileName)
{    //printf("++++++++++CALL TO PRINT THE LIST__________+++++++++++++++++++\n");
    if(isDataLstEmpty(dataHdrInfo))
        return;
    fflush(stdout);
    FILE *fp=fopen(strFileName,"w");
    //printf("++++++++++CALL TO PRINTing THE representative  LIST +++++++++++++++++++\n");
    int iCnt = 0,j;
    //DataNd dataTemp = dataHdrInfo->dataFirstNd;
    for(j=0;j<dataHdrInfo->uiCnt;j++)
    {
      //  fprintf(fp,"%d %d %d\n",dataHdrInfo->dataClstElem[j].iNum,dataHdrInfo->dataClstElem[j].core_tag,dataHdrInfo->dataClstElem[j].dataCluster);
        
        /*for(iCnt = 0; iCnt < DIMENSION; iCnt++ )
            printf("%lf ", dataHdrInfo->dataClstElem[j].iData[iCnt] );
            */
        

    }
    fclose(fp);

  return;
}


void printData(Data dataPoint)
{
    printf("RecNum %d\t",dataPoint->iNum);
    int iCnt=0;
        for(iCnt = 0; iCnt < DIMENSION; iCnt++)
            printf("%lf ", dataPoint->iData[iCnt] );

    return;
}

Boolean isDataLstEmpty(DataHdr dataHdrInfo)
{
    return(dataHdrInfo == NULL || dataHdrInfo->uiCnt == 0) ? TRUE : FALSE;
}

void freeDataList(DataHdr dataList1)
{
    int i;
    for(i=0;i<dataList1->uiCnt;i++)
        free((dataList1->dataClstElem+i)->iData);
    free(dataList1->dataClstElem);
    free(dataList1);

    return;
}

void writeDataListToFile(DataHdr dataList1, char * fileName)
{

    FILE *filP = fopen(fileName, "w");

    int i =0; int j=0;

    fprintf(filP, "%d\n", dataList1->uiCnt);
    fprintf(filP, "%d\n", DIMENSION);

    for(i=0;i<dataList1->uiCnt;i++)
    {
        for (j = 0; j < DIMENSION; j++)
        {
            fprintf(filP, "%lf", dataList1->dataClstElem[i].iData[j]);
            if(j!=(DIMENSION-1))
            {
                fprintf(filP," ");
            }
        }

        fprintf(filP,"\n");
    }

    fclose(filP);

    return;

}
