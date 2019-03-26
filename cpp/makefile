FLAG=-c -g
CC=g++

output: Data.o Grid.o GList.o RList.o BCell.o RTree.o GridRTree.o main.o 
	$(CC) -o output Data.o Grid.o GList.o RList.o BCell.o RTree.o GridRTree.o main.o -lm  

main.o: main.cpp
	$(CC) $(FLAG) main.cpp -lm

Grid.o: Grid.cpp
	$(CC) $(FLAG) Grid.cpp 

Data.o: Data.cpp
	$(CC) $(FLAG) Data.cpp -lm
	
GridRTree.o: GridRTree.cpp
	$(CC) $(FLAG) GridRTree.cpp -lm

GList.o: GList.cpp
	$(CC) $(FLAG) GList.cpp -lm 

RList.o: RList.cpp
	$(CC) $(FLAG) RList.cpp -lm 

BCell.o: BCell.cpp
	$(CC) $(FLAG) BCell.cpp -lm 

RTree.o: RTree.cpp
	$(CC) $(FLAG) RTree.cpp -lm

clean:
	rm *.o
	rm output
