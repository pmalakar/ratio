#ifndef __personality__
#define __personality__

#include <mpix.h>
#include <stdint.h>

//Variables
/* 
 * Size of midplane on BGQ
 */
extern const int MidplaneSize; 

/*
 * Number of midplanes : Number of processes / MidplaneSize
 */
extern int numMidplanes, midplane;

/* 
 * Structure containing BGQ parameters 
 */
extern MPIX_Hardware_t hw; 

/*
 * Ranks per node, Core ID (0...15) 
 */ 
extern int ppn, coreID, nodeID;

/*
 * Size of each dimension 
 */ 
extern int dimSize[MPIX_TORUS_MAX_DIMS];	// MPIX_TORUS_MAX_DIMS=torus dimension size

/*
 * Whether dimension is torus or not	//TODO bool 
 */ 
extern int isTorus[MPIX_TORUS_MAX_DIMS];	// torus wrap = 0/1 

/*
 *  2D array for each rank, contains bridge node info
 *  [0] - Rank of the bridge node, [1] - Distance from the IO node 
 */ 
extern int bridgeNodeInfo[2];						

/*
 *  Routing order : varies based on the partition, node count etc. 
 */
extern int *routingOrder;

/*
 * 	Maximum heap available per rank 
 */
extern uint64_t heapAvail;

//Functions
void getPersonality(int, int);

#endif
