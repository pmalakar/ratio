#include <stdio.h>
#include <assert.h>
#include <mpi.h>
#include <mpix.h>

#include "personality.h"
#include "neighbour.h"

//Variables
int **neighbourRanks;	//[MidplaneSize][10];
int NUM_NEIGHBOURS=10;

//called before calling findNeighbours for the first time
void initNeighbours(int size) {
	neighbourRanks = new int * [size]; // [MidplaneSize*ppn];
	for (int i=0 ; i<size ; i++)
		neighbourRanks[i] = new int[NUM_NEIGHBOURS];
}

/*
	* change each coordinate of the current rank to get ten neighbours
	 
	//neighbour (2D array) of (a, b, c, d, e)
	//neighbour [0]  a+ b c d e t
	//neighbour [1]  a- b c d e t
	//neighbour [2]  a b+ c d e t
	// ...
*/
int findNeighbours (int myrank, int currentRank) {

	// * * * * * Variable declaration * * * * * //
 
	int i, rank, k=0, count = -1, num, dim, index;
	int currentCoords[MPIX_TORUS_MAX_DIMS+1];
	int neighbour[10][6];	//max 10 neighbours (10 links)

	assert(neighbourRanks);

	//Get coordinates of current node	
	MPIX_Rank2torus (currentRank, currentCoords);

	//Initialize ranks of neighbours
	for (i=0; i<10; i++)
		neighbourRanks[currentRank][i] = -1;

	//Initialize neighbour coordinates to current rank coordinates
	for (num=0; num<10; num ++) 
	 for (dim=0; dim<=MPIX_TORUS_MAX_DIMS; dim++) { 
		neighbour[num][dim] = currentCoords[dim];
	 }
	
	//Deduce neighbour coordinates and ranks, traverse plus and minus directions of each dimension (A, B, C, D, E)
	for (dim=0; dim < MPIX_TORUS_MAX_DIMS ; dim++) { 

		// plus direction

		index = 2*dim;
#ifdef DEBUG
		printf("1. %d mycoord[%d]=%d neighbour #%d nghbrdim now =%d\n", currentRank, dim, currentCoords[dim], index, neighbour[index][dim]);
#endif
		/*if (currentCoords[dim] < hw.Size[dim]-1)
		 neighbour[index][dim] = currentCoords[dim] + 1;
		else {
			if (hw.isTorus[dim] == 1) neighbour[index][dim] = 0;
			else neighbour[index][0] = -1, neighbour[index][dim] =-1;	//mesh, no neighbour in this direction
		}*/

		//TODO fix...
		//if (hw.isTorus[dim] == 1) neighbour[index+1][dim] = (currentCoords[dim] - 1 + hw.Size[dim]) % hw.Size[dim];		//torus
		if (hw.isTorus[dim] == 1 || (currentCoords[dim] < hw.Size[dim]-1))		//either torus or not the last row or col in mesh	
			neighbour[index][dim] = (currentCoords[dim] + 1) % hw.Size[dim];	// torus
		else
			neighbour[index][0] = -1, neighbour[index][dim] =-1;				// mesh, no neighbour in this direction
		
#ifdef DEBUG
		printf("2. %d dim=%d nghbr %d nghbrdim new =%d\n", currentRank, currentCoords[dim], index, neighbour[index][dim]);
		printf("2. for %d nghbrdim new coord = %d %d %d %d %d\n", currentRank, neighbour[index][0],neighbour[index][1],neighbour[index][2],neighbour[index][3],neighbour[index][4]);
#endif

		//Get the rank of neighbour[index]
		if(neighbour[index][0] != -1 && neighbour[index][dim] != -1) { 
		 
		  MPIX_Torus2rank(neighbour[index], &rank); 
		  neighbourRanks[currentRank][++count] = rank;

#ifdef DEBUG
		  printf("%d: Set neighbour[%d] for %d as %d i.e. %d (%d %d %d %d %d)\n", \
			myrank, count, currentRank, neighbourRanks[currentRank][count], rank, \
			neighbour[index][0], neighbour[index][1], neighbour[index][2], neighbour[index][3], neighbour[index][4]);
#endif
		}
		else
			neighbourRanks[currentRank][++count] = -1;

		// minus direction
		// if dim size is 2, there is only 1 direction
		if (dimSize[dim] == 2)
			continue;
		
		//index+1
		//k=0;
#ifdef DEBUG
		printf("%d: 3. %d mycoord[%d]=%d neighbour #%d nghbrdim now =%d\n", myrank, currentRank, dim, currentCoords[dim], index+1, neighbour[index+1][dim]);
#endif
		/*if (currentCoords[dim] > 0)
			neighbour[index+1][dim] = currentCoords[dim] - 1;
		else {
			if (hw.isTorus[dim] == 1) neighbour[index+1][dim] = hw.Size[dim]-1;
			else neighbour[index+1][0] = -1, neighbour[index+1][dim] =-1;	//mesh
		}*/
		
		//testing plus direction next ... for plus neighbours being found first from the reverse 
		//TODO fix...
		//if (hw.isTorus[dim] == 1)		neighbour[index][dim] = (currentCoords[dim] + 1) % hw.Size[dim];	// torus
		if (hw.isTorus[dim] == 1 || currentCoords[dim] > 0) 
			neighbour[index+1][dim] = (currentCoords[dim] - 1 + hw.Size[dim]) % hw.Size[dim];		//torus
		else 
			neighbour[index+1][0] = -1, neighbour[index+1][dim] =-1;							//mesh

#ifdef DEBUG
		printf("%d: 4. %d dim=%d neighbour %d nghbrdim new =%d\n", myrank, currentRank, currentCoords[dim], index+1, neighbour[index+1][dim]);
		printf("%d: 4. for %d nghbrdim new coord = %d %d %d %d %d\n", myrank, currentRank, neighbour[index+1][0],neighbour[index+1][1],neighbour[index+1][2],neighbour[index+1][3],neighbour[index+1][4]);
#endif

		//get the rank of neighbour[index]
		if(neighbour[index+1][dim] != -1) { 

		  MPIX_Torus2rank(neighbour[index+1], &rank);

		  //check in neighbourRanks, if rank is not present then add
		  for (k=0; k<=count; k++) {
#ifdef DEBUG
			printf("%d: compare for %d at %d ---- %d %d\n", myrank, currentRank, k, neighbourRanks[currentRank][k], rank);
#endif
			if(rank == neighbourRanks[currentRank][k]) { 
#ifdef DEBUG
				printf("%d: Rank %d at index+1 %d was already added for %d\n", myrank, rank, index+1, currentRank); break; 
#endif
			}
		  }

		  if (k>count) {	//no matching rank found, add the rank 
				neighbourRanks[currentRank][++count] = rank;
#ifdef DEBUG
				printf("%d: Set neighbour[%d] for %d as %d i.e. %d (%d %d %d %d %d)\n", myrank, count, currentRank, neighbourRanks[currentRank][count], rank, neighbour[index+1][0], neighbour[index+1][1], neighbour[index+1][2], neighbour[index+1][3], neighbour[index+1][4]);
#endif
			}
		}
		else
			neighbourRanks[currentRank][++count] = -1;
	}

#ifdef DEBUG
	for (num = 0; num<10; num ++) {
	  if (neighbourRanks[currentRank][num] != -1) {
			printf ("%d: Neighbour %d of %d (%d %d %d %d %d) at %d => %d %d %d %d %d => %d\n", myrank, num, currentRank, currentCoords[0], currentCoords[1], currentCoords[2], currentCoords[3], currentCoords[4], bridgeNodeInfo[0], neighbour[num][0], neighbour[num][1], neighbour[num][2], neighbour[num][3], neighbour[num][4], neighbourRanks[currentRank][num]);
	  }	
	}
	
	for (num = 0; num<10; num ++) 
	  if (neighbourRanks[currentRank][num] != -1)  
			printf (" %d -- %d;\n", currentRank, neighbourRanks[currentRank][num]);
#endif
	
	return 0;
}


