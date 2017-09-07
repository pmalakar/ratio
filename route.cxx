#include "route.h" 

/*
 * Read the DCR registers to get the routing order 
 * The routing order depends on the partition in which the job runs
 * [Thanks to Phil Heidelberger]
 *
 * getRoutingOrder() populates parameter ro with the routing order, viz.  
 * A=0 B=1 C=2 D=3 E=4
 *
 * Example, for 512 nodes on Mira, ro = {0, 1, 2, 3, 4}
 *
 */

int getRoutingOrder (int *ro) {

	uint64_t dcr_det_order = DCRReadUser(ND_500_DCR(CTRL_DET_ORDER));

	int A = ND_500_DCR__CTRL_DET_ORDER__MASK0_get(dcr_det_order);
	int B = ND_500_DCR__CTRL_DET_ORDER__MASK1_get(dcr_det_order);
	int C = ND_500_DCR__CTRL_DET_ORDER__MASK2_get(dcr_det_order);
	int D = ND_500_DCR__CTRL_DET_ORDER__MASK3_get(dcr_det_order);
	int E = ND_500_DCR__CTRL_DET_ORDER__MASK4_get(dcr_det_order);

	int torus[5] = {A, B, C, D, E};

	int index = 0;
  for (int i=0; i<5 ; i++) { 
	 	if (torus[i] == 1) 					index = 4;
	 	else if (torus[i] == 2)			index = 3;
	 	else if (torus[i] == 4)			index = 2;
	 	else if (torus[i] == 8)			index = 1;
	 	else if (torus[i] == 16)		index = 0;
    ro[i] = index;
	}

	return 0;
}

/*
int main (int argc, char *argv[]) {

  int myrank;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

  if (myrank == 0) {
  	int *ro = new int[5];
		getRouting(ro);
 		printf("%d: %d %d %d %d %d\n", myrank, ro[0], ro[1], ro[2], ro[3], ro[4]);
  	delete[] ro;
	}

	MPI_Finalize ();
  
}
*/

