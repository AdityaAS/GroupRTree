#include "Def.h"
#include "GridRTree.h"
#include "Grid.h"
#include "BCell.h"
#include "Data.h"
#include <ctime>
#include <iostream>
#include <vector>
#include <queue>
#include <deque>
#include <cstdlib>

using namespace std;

#define CORE_POINT 1
#define NOT_CORE_POINT 0

#define SUCCESS 0
#define FAILURE -3
#define NOISE -1
#define NOGROUP -5

struct Cluster
{
	int cluster_id;
	int core_pts;
	vector<int> point_ids;
};

int countGR, countDB;
int *arrDB, *arrGR, *cluster1, *cluster2;
FILE *joining1 , *joining2;

double EPS;
double EPS1;
int MINPOINTS;
int UNDEFINED;
int DIMENSION;    //value for the number of dimensions
int noOfPoints;
int GMINENTRIES;  
int GMAXENTRIES;
int GAUXMINENTRIES;
int GAUXMAXENTRIES; 
int RMINENTRIES;        //Minimum entries for node in RTree
int RMAXENTRIES;       //Maximum entries for node in RTree
double CELLSIZE;
int TIMESTAMP;
int CHOICE;
int CLUSTERID;
int IMMCOUNT;
int EPSCOUNT;
int noise;
int BCELLID;
int GROUPID;
long long int DATAID;
double * MINGRIDSIZE;
double * MAXGRIDSIZE;
double *MINGRIDSIZEglobal, *MAXGRIDSIZEglobal;
double **MinGridBound, **MaxGridBound;
vector<Group> groupList;
int totalPoints;

FILE* clusterOuputFile;

void dbscan(DataHdr dataList, double epsilon, int min_pts);
int expand(int index, int cluster_id, DataHdr dataList, int num_points, double epsilon, int min_pts);

void dbscan(DataHdr dataList, double epsilon, int min_pts)
{
    int cluster_id = 0;
    for (int i = 0; i < totalPoints; ++i) 
    {
    	Data points = dataList->dataClstElem + i;
        if (points->ClusterID == UNDEFINED) 
        {
            if (expand(i, cluster_id, dataList, totalPoints, epsilon, min_pts) == CORE_POINT)
            {
                cluster_id+=1;
            }
        }
    }

    Cluster** clusters = (Cluster**)malloc(sizeof(Cluster*)*cluster_id);
    vector<int> clusterpoints;

    for(int i=0;i<cluster_id;i++)
    {
        clusters[i] = new Cluster;
        clusters[i]->cluster_id = i;
        clusters[i]->core_pts = 0;
        clusters[i]->point_ids.clear();
        clusterpoints.push_back(0);
    }

    noise = 0;
    for(int i=0;i<totalPoints;i++)
    {
    	Data datapoint = dataList->dataClstElem + i;

        if(datapoint->ClusterID == UNDEFINED)
        {
            cout << "Cluster Unassigned Error" << endl;
            return;
        }

        else if(datapoint->ClusterID == NOISE)
        {
            noise++;
        }

        else
        {
            clusters[datapoint->ClusterID]->point_ids.push_back(i);

            if(datapoint->core_tag == TRUE) clusters[datapoint->ClusterID]->core_pts++;
            
            clusterpoints[datapoint->ClusterID]++;
        }
    }
    
    // cout << "Noise: " << noise << endl;
    // int border = 0;
    int core = 0;
    for(int i=0;i<cluster_id;i++)
    {
        // border = border + (clusters[i]->point_ids.size() - clusters[i]->core_pts);
        core = core + clusters[i]->core_pts;
        // fprintf(stderr, "Cluster: %d \tTotal: %d \t Core: %d\n", i, clusters[i]->point_ids.size(), clusters[i]->core_pts);
        // cout << "Cluster: " << i << "\tTotal: " << clusters[i]->point_ids.size() << "\tCore: " << clusters[i]->core_pts << endl;
    }
    // fprintf(stderr, "", core);
    fprintf(stderr, "Clusters: %d\tCore:%d\t Noise:%d\n", cluster_id, core, noise);
    // cout << "No of Clusters: " << cluster_id << endl;
    // fprintf(stderr, "Noise:\t %d\n", noise);
    clusterpoints.clear();

    for(int i=0;i<cluster_id;i++)
    {
    	clusters[i]->point_ids.clear();
    	free(clusters[i]);
    }
    free(clusters);
}

int expand(int index, int cluster_id, DataHdr dataList, int num_points, double epsilon, int min_pts)
{
    int return_value = NOT_CORE_POINT;
    vector<int> seeds = findNeighbours(dataList, index, epsilon);
    Data currPoint = dataList->dataClstElem + index;
    Group g = groupList[currPoint->group_id];

    if (seeds.size() < min_pts)
    {
        currPoint->ClusterID = NOISE;
    }
    else
    {
    	currPoint->core_tag = TRUE;
        currPoint->ClusterID = cluster_id;
        queue<int, deque<int> > q (deque<int> (seeds.begin(),seeds.end())) ;

        for(int i=0;i<seeds.size();i++)
            (dataList->dataClstElem + seeds[i])->ClusterID = cluster_id;

        while(!q.empty())
        {
            int id = q.front();
            q.pop();
            vector<int> spread = findNeighbours(dataList, id, epsilon);
            
            if(spread.size() >= min_pts) 
            {
            	Data point = dataList->dataClstElem + id;
            	point->core_tag = TRUE;
                // points[id]->type = CORE;
                for(int i=0;i<spread.size();i++)
                {
                    if ((dataList->dataClstElem + spread[i])->ClusterID == NOISE || (dataList->dataClstElem + spread[i])->ClusterID == UNDEFINED) 
                    {
                        if ((dataList->dataClstElem + spread[i])->ClusterID == UNDEFINED)
                        {
                            q.push((dataList->dataClstElem + spread[i])->id);
                        }
                        (dataList->dataClstElem + spread[i])->ClusterID = cluster_id;
                    }
                }
            }
            spread.clear();
        }

        return_value = CORE_POINT;
    }
    seeds.clear();
    return return_value;
}

int main(int argc, char **argv){

    clock_t global_t = clock();
    clock_t local_t;
    int myrank, numprocs, i, j;
    
    printf("%d\n", argc);

    if(argc != 7)
    {
        fprintf(stderr, "Usage: ./<program> <inputFile> <m of RTree> <M of RTree> <Epsilon> <MinPts> <ClusterOutput>\n");
        // cout << "Usage: ./<program> <inputFile> <m of RTree> <M of RTree> <Epsilon> <MinPts> <ClusterOutput>" << endl;
        return 0;
    }
	struct stat st;
	int fileExists = stat(argv[1], &st);
    fprintf(stderr, "Input: %s\tEps: %lf\t Minpts: %d\n", argv[1], strtod(argv[4],NULL), atoi(argv[5]));
	if (fileExists != 0){
		#ifdef PRINT
			if(myrank == 0)printf("Input file doesn't exist.\n");
		#endif
		exit(-1);
	}

    GMINENTRIES = atoi(argv[2]);
    GMAXENTRIES = atoi(argv[3]);
    EPS = strtod(argv[4],NULL);
    MINPOINTS = atoi(argv[5]);

    GAUXMINENTRIES = GMINENTRIES;
    GAUXMAXENTRIES = GMAXENTRIES;

    RMINENTRIES = GMINENTRIES;
    RMAXENTRIES = GMAXENTRIES;

	CLUSTERID=0;
	BCELLID=0;
	GROUPID=0;

	UNDEFINED = 100000000;

	DataHdr dataList;
	GHdrNd GRTree;

	local_t = clock();
	dataList = readData(argv[1]);
	local_t = clock() - local_t;
    // fprintf(stderr, "Time for Reading Input\t %Lf\n", ((long double)local_t)/CLOCKS_PER_SEC);
    local_t = clock();
	
	CELLSIZE = EPS/ (2*sqrt(DIMENSION));
	GroupListHd groupListhd = initGroupListHd();

	GRTree = populateGridRTree(dataList, groupListhd, GMINENTRIES, GMAXENTRIES);
	
	constructAuxRTrees(dataList, groupListhd);
	local_t = clock() - local_t;
    // fprintf(stderr, "Time for Group R Tree Construction\t %Lf\n", ((long double)local_t)/CLOCKS_PER_SEC);
    local_t = clock();
	totalPoints=noOfPoints;

	findReachableGroupsRTree(GRTree, groupList);
	local_t = clock() - local_t;
    long double reachabletime;
    reachabletime = (long double)local_t/CLOCKS_PER_SEC;
    // fprintf(stderr, "Time for Reachable Groups Calculation\t %Lf\n", ((long double)local_t)/CLOCKS_PER_SEC);
    local_t = clock();
	dbscan(dataList, EPS, MINPOINTS);
	// freeDataList(dataList);

 // 	freeGRTree(GRTree);

 //    freeGroupListHd(groupListhd);

    groupList.clear();
	local_t = clock() - local_t;
    // fprintf(stderr, "Time for DBSCAN\t %Lf\n", ((long double)local_t)/CLOCKS_PER_SEC);

    global_t = clock() - global_t;
    // fprintf(stderr, "Total Runtime\t %Lf\n", ((long double)local_t)/CLOCKS_PER_SEC);

    return 0;
}

	
