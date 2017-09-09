#ifndef __common__
#define __common__

#include <iostream>
#include <iomanip>

#define MAXBUF (1024*32)

using namespace std;

class dataBlock {

 private:	
	double *alpha;
	int numElem;

 public:
	dataBlock() {}
	dataBlock(int count) {
		numElem = count;
	}

	void allocElement () {
		try {
			alpha = new double[numElem];
			for (int i=0; i<numElem ; i++) {
				alpha[i] = rand() % 100;
			}
		}
		catch (bad_alloc& ba) {
				cerr << "Bad allocation for alpha\n" << ba.what() << endl;
		}
	}

	void freeElement () {
		delete [] alpha;
	}

	double *getAlphaBuffer() {
		return &alpha[0];
	}

};


inline void prnerror (int error_code, char *string)
{
	
	char error_string[256];
	int length_of_error_string;	
	MPI_Error_string(error_code, error_string, &length_of_error_string);
	fprintf(stderr, "%3d: %s in %s\n", error_code, error_string, string);
	MPI_Finalize();
	exit(-1);
}

#endif
