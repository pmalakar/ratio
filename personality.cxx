/*
 *  Created by Preeti Malakar
 *  Argonne National Laboratory
 *
 *  For BG/Q: Determine routing order, bridge node information etc.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <mpi.h>
#include <mpix.h>

#ifdef BGQ
#include <spi/include/kernel/memory.h>
#endif

#include "personality.h"
#include "route.h"

/* 
 * Size of midplane on BGQ
 */
const int MidplaneSize = 512; 

/*
 * Number of midplanes : Number of processes / MidplaneSize
 */
int numMidplanes, midplane;

/* 
 * Structure containing BGQ parameters 
 */
MPIX_Hardware_t hw; 

/*
 * Ranks per node, Core ID (0...15) 
 */ 
int ppn, coreID, nodeID;

/*
 * Size of each dimension 
 */ 
int dimSize[MPIX_TORUS_MAX_DIMS];	// torus dimension size

/*
 * Whether dimension is torus or not	//TODO bool 
 */ 
int isTorus[MPIX_TORUS_MAX_DIMS];	// torus wrap = 0/1 

/*
 *  2D array for each rank, contains bridge node info
 *  [0] - Rank of the bridge node, [1] - Distance from the bridge node 
 */ 
int bridgeNodeInfo[2];						

/*
 *  Routing order : varies based on the partition, node count etc. 
 */
int *routingOrder;

/*
 *	Maximum heap available per rank 
 */
uint64_t heapAvail;

inline int min (int x, int y) {
	return x<y? x: y;
}


/*
 * Initialize system parameters
 * - get ranks per node, nodeID, coreID, torus wrap=0/1?
 * - get routing order
 */

void initSystemParameters(int rank) {

	int i;

	routingOrder = new int[MPIX_TORUS_MAX_DIMS];
	getRoutingOrder(routingOrder);

	MPIX_Hardware(&hw);
	ppn = hw.ppn;
	coreID = hw.coreID;	
	nodeID = rank/ppn;

	for (i=0; i<MPIX_TORUS_MAX_DIMS ; i++) {
		isTorus[i] = hw.isTorus[i];
		dimSize[i] = hw.Size[i];
	}

	//Rank of bridge node and distance to IO node
	bridgeNodeInfo[0] = MPIX_IO_link_id (); 
	bridgeNodeInfo[1] = MPIX_IO_distance (); 

#ifdef STATS
#ifdef BGQ
	Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapAvail); 
#endif
#endif

}

/*
 * if destRank is -1, this function call was just made to populate bridgenodeinfo
 *
 * if not, then this function was used to print the path from rank to destRank (in order to visualize and debug)
 *
 */

void getPersonality (int rank, int destRank) {

	//local variables
	int i;
	int unitHop = 1;

	int flag=0;

	if (destRank == -1) {
		initSystemParameters(rank);
		destRank = bridgeNodeInfo[0];
		flag = 1;
	}

#ifdef DEBUG

	printf("Torus dimensions = (%u,%u,%u,%u,%u) Routing order = (%d,%d,%d,%d,%d)\n", hw.Size[0], hw.Size[1], hw.Size[2], hw.Size[3], hw.Size[4], routingOrder[0], routingOrder[1], routingOrder[2], routingOrder[3], routingOrder[4]);

	if (rank == 0)
		printf("Torus wraps? %u,%u,%u,%u,%u\n", hw.isTorus[0], hw.isTorus[1], hw.isTorus[2], hw.isTorus[3], hw.isTorus[4]);

	if (flag == 1)
		printf("Rank: %d Node: %d Torus coords = (%u,%u,%u,%u,%u) distance to ION: %d link ID: %d\n", rank, rank/ppn, hw.Coords[0], hw.Coords[1], hw.Coords[2], hw.Coords[3], hw.Coords[4], bridgeNodeInfo[1], destRank);

#endif

	if (flag == 1)
		return;

	/*
	 * Find coordinates of bridge node
	 */
	int destCoords[6];
	MPIX_Rank2torus (destRank, destCoords);

	/*
	 * Initialize intermediate nodes in original path to the bridge node
	 */
	int intmdtCoords[6];
	for (int dim=0; dim < MPIX_TORUS_MAX_DIMS; dim++) 
		intmdtCoords[dim] = hw.Coords[dim];

	intmdtCoords[MPIX_TORUS_MAX_DIMS] = 0;

	int hopnum = 0;
	int hopDiff, intmdt_rank, child, parent;
	child = rank;
	for (int dim=0; dim<MPIX_TORUS_MAX_DIMS; dim++) {

		int dimID = routingOrder[dim];
		hopDiff = abs(destCoords[dimID] - hw.Coords[dimID]);

		//	if (hw.isTorus[dimID] == 1 && (hopDiff*2 > hw.Size[dimID])) 
		//			hopDiff = hw.Size[dimID] - hopDiff ;

		if (hw.isTorus[dimID] == 1) 
			hopDiff = min (hopDiff, hw.Size[dimID] - hopDiff) ;

#ifdef DEBUG
		if (flag == 0)
			printf("%d to %d difference in dim %d = %d\n", rank, destRank, dimID, hopDiff);
#endif

		for(int diff=0; diff<hopDiff ;diff++) {
			if (hw.isTorus[dimID] == 0) {
				if(destCoords[dimID] < hw.Coords[dimID]) intmdtCoords[dimID] -= unitHop;  
				else intmdtCoords[dimID] += unitHop;
			}
			else {		// torus
				if (abs(destCoords[dimID] - hw.Coords[dimID])*2 > hw.Size[dimID]) {
					//	printf("check > %d bridgecoords[%d]=%d hw.Coords[%d]=%d hw.Size[%d]=%d\n", rank, dimID, destCoords[dimID], dimID, hw.Coords[dimID], dimID, hw.Size[dimID]);

					if (destCoords[dimID] > hw.Coords[dimID]) 
						intmdtCoords[dimID] = ((intmdtCoords[dimID] - unitHop) + hw.Size[dimID]) % hw.Size[dimID];  
					else 
						intmdtCoords[dimID] = (intmdtCoords[dimID] + unitHop) % hw.Size[dimID];
				}
				else if (abs(destCoords[dimID] - hw.Coords[dimID])*2 < hw.Size[dimID]) {
					//	printf("check < %d bridgecoords[%d]=%d hw.Coords[%d]=%d hw.Size[%d]=%d\n", rank, dimID, destCoords[dimID], dimID, hw.Coords[dimID], dimID, hw.Size[dimID]);
					if (destCoords[dimID] < hw.Coords[dimID]) 
						intmdtCoords[dimID] = ((intmdtCoords[dimID] - unitHop) + hw.Size[dimID]) % hw.Size[dimID];  
					else 
						intmdtCoords[dimID] = (intmdtCoords[dimID] + unitHop) % hw.Size[dimID];
				}
				else {
					//if source coord is even, plus direction
					if (hw.Coords[dimID]%2 == 0)	// see phil's email: Aug 22, 2014
						intmdtCoords[dimID] = (intmdtCoords[dimID] + unitHop) % hw.Size[dimID];			//even source coord: traverse in plus direction  
					else 
						intmdtCoords[dimID] = ((intmdtCoords[dimID] - unitHop) + hw.Size[dimID]) % hw.Size[dimID];  
				}
			}

			++hopnum;

			//get the intermediate node rank
			MPIX_Torus2rank (intmdtCoords, &intmdt_rank);
			parent = intmdt_rank;

			if (flag == 0) {			
				printf ("Enroute %d (%d %d %d %d %d) to %d (%d %d %d %d %d) Hop %d: in dimension %d Child %d to Parent %d (%d %d %d %d %d)\n", \
						rank, hw.Coords[0], hw.Coords[1], hw.Coords[2], hw.Coords[3], hw.Coords[4], \
						destRank, destCoords[0], destCoords[1], destCoords[2], destCoords[3], destCoords[4], \
						hopnum, dimID, child, intmdt_rank, intmdtCoords[0], intmdtCoords[1], intmdtCoords[2], intmdtCoords[3], intmdtCoords[4]);

#ifdef DEBUG
				printf ("Route %d to %d Hop %d\n", rank, intmdt_rank, hopnum);
#endif
				printf ("%d->%d;\n", child, parent);

			}

			child = parent;
		}
	}   

#ifdef DEBUG
	//Check if everyone was routed to their bridge nodes?
	if (parent != destRank && bridgeNodeInfo[1] > 1) 
		printf("Rank %d lost, did not reach %d, instead captivated by %d :D\n", rank, destRank, parent);
	if (hopnum != bridgeNodeInfo[1]-1)  
		printf("Rank %d differs in number hops %d!=%d !!!!\n", rank, bridgeNodeInfo[1], hopnum);
	printf ("%d reaches %d in %d hops\n", rank, destRank, hopnum);
#endif

	return;

}


