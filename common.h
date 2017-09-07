#ifndef __common__
#define __common__

#define MAXBUF (1024*32)

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

