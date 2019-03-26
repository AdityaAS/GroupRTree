#ifndef __DATA_H
#define __DATA_H

#include "Def.h"

void insertDataLstElem(DataHdr dataHdrInfo, DataPoint iData);
void insertDataLstElem1(DataHdr dataHdrInfo, DataPoint iData, int index, int serial);
DataHdr initDataHdr(int size_data);
Data initData(DataPoint iData);
DataHdr readData(char *strFileName);
int readPoints(char *strFileName);
DataHdr readSampledData(char *strFileName);
void printDataLst(DataHdr dataHdrInfo);
void printDataLstFile(DataHdr dataHdrInfo,char * strFileName);
void printData(Data dataPoint);
Boolean isDataLstEmpty(DataHdr dataHdrInfo);
void freeDataList(DataHdr dataList1);
void writeDataListToFile(DataHdr dataList1, char * fileName);
void deleteData(DataHdr dataList);
void freeRNbHdr(RNbHdr r);

#endif