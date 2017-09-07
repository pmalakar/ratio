#ifndef __neighbour__
#define __neighbour__

//neighbour list of all nodes
extern int **neighbourRanks;// = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
//extern int neighbourRanks[][10];

void initNeighbours (int); 
int findNeighbours (int, int); 

#endif
