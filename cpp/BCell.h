#ifndef BCELL_H
#define BCELL_H

#include "Def.h"

GroupListHd initGroupListHd();
Group initGroup(Data datapoint);
CellDataHd initCellDataHd();
BCell initBCell(Region rectTemp);
CellData initCellData(Data dataClstElem);
BCellListHd initBCellListHd();
void freeCellsList(BCellListHd cellsListNbh);
void freeAllBCells(BCellListHd cellsList);
void freeBCell(BCell bCellElem);
void freeCellDataList(CellDataHd cellDataList);
void freeGroupListHd(GroupListHd groupList);

#endif
