#ifndef __iotree__
#define __iotree__

#define NONBLOCKING 0
#define prnl printf("\n")
#define __STDC_FORMAT_MACROS 1 // In C++ the macros are not automatically defined just by including the file. 

#define MAXBUF (1024*32)
#define xMicroSec 1000000

using namespace std;

// Variables
extern int myrank, commsize, mode, fileSize;

class Node {

	private:
		int nodeId, rootId, lastIndex;
		int depth;				//depth of the node in the tree
		int childIndex[9];
		Node *childIndexPtr[9];
		Node *parent;

	public:

		Node *leftSibling, *rightSibling;	//tree
		Node *prev, *next;					//queue
		Node() { lastIndex = -1;}
		~Node() {}

		Node (int currentNode) {
			lastIndex = -1;
			nodeId = currentNode;			
			parent = NULL;
			depth=0;
			rootId = currentNode;
			//printf("depth %d\n", depth);
		}

		Node* addChild (int id, int rootid) {
			Node *node = new Node(id);
			node->parent = this;
			node->depth = depth + 1;
			node->rootId = rootid;
			childIndexPtr[++lastIndex] = node;
			return node;
		}		

		Node* getChild (int i) {
			if (i>lastIndex) return NULL;
			return childIndexPtr[i];
		}
	
		int getChildId (int i) {
			if (i>lastIndex) return -1;
			return childIndexPtr[i]->getNodeId();
		}

		int getNodeId () 			{ return nodeId; }
		int getRootId () 			{ return rootId; }
		int getNumChildren () { return lastIndex+1; }
		int getDepth() 				{ return depth; }
		Node *getParent () 		{ return parent; }

		void printChildrenInfo (int rank) {
			int i;
			printf("%d: num children of %d = %d\n", rank, nodeId, lastIndex);
			for (i=0; i<=lastIndex; i++)
				if (rank == nodeId) printf("First level of (%d): %d[%d]\n ", nodeId, childIndexPtr[i]->getNodeId(), i);
		}

};


//data used for computation, analysis, IO etc..

class dataBlock {

 private:	
	double *alpha, *beta, *gamma;
//    double U[1024][1024], V[1024][1024], W[1024][1024];
	double **U, **V, **W;
	int numElem;
	int latsize, lonsize;

 public:
	dataBlock() {}
	dataBlock(int count) {

		numElem = count;

	}

	void allocElement (int type) {
		
		if (type == 1) {
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
	}

	void freeElement (int type) {

		if (type == 1) delete [] alpha;
	}

	double *getAlphaBuffer() {
		return &alpha[0];
	}

};

// Function Declarations

int writeFile (dataBlock *, int, int);
void coalesceData (dataBlock *, int);
void getData(int);
int findNeighbours (int);

inline void prnerror (int error_code, char *string)
{
	
	char error_string[256];
	int length_of_error_string;	
	MPI_Error_string(error_code, error_string, &length_of_error_string);
	fprintf(stderr, "%3d: %s in %s\n", error_code, error_string, string);
	MPI_Finalize();
	exit(-1);
}

inline int min (int a, int b) {
	if (a<=b) return a;
	else return b;
}

#endif
