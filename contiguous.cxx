#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include "common.h"
#include "algorithm.h"

#define oneKB 1024
#define NAME_LENGTH 256

int rank, size;
int count, mode;

double tION, tFS;

MPI_File fileHandle;
MPI_Status status;
MPI_Request request;

char *fileNameION = "/dev/null";
char *fileNameFS = "dummyFile";
char *fileNameFSBN = "dummyFileBN";
char *fileNameFSCO = "dummyFileCO";
char testFileName[NAME_LENGTH];

int SKIP = 1;
int MAXTIMES = 5;
int MAXTIMESD = 7;

int writeFile(dataBlock *datum, int count) 
{
  int result, nbytes;
  int blocking = 0;

	MPI_Request request;

	if (blocking == 1) {
		result = MPI_File_write_at (fileHandle, (MPI_Offset)rank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
		if (result != MPI_SUCCESS) 
			prnerror (result, "all MPI_File_write_at Error:");
	}
	else {
		MPIO_Request req_iw;
		MPI_Status st_iw;
		result = MPI_File_iwrite_at (fileHandle, (MPI_Offset)rank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &req_iw);
		MPI_Wait(&req_iw, &st_iw);
	}

	MPI_Get_elements( &status, MPI_CHAR, &nbytes );

	return nbytes;
}

void file_write(dataBlock *datum) {

  int totalBytes = 0;
  double tIOStart, tIOEnd;

#ifdef BGQ
	/*
	 * * * * * * * * * Independent MPI-IO to IO nodes from all compute nodes - shared file * * * * * * * *
   */
  MPI_File_open (MPI_COMM_WORLD, fileNameION, mode, MPI_INFO_NULL, &fileHandle);
  for (int i=1; i<=SKIP; i++)
   totalBytes += writeFile(datum, count);
  tIOStart = MPI_Wtime();
  for (int i=1; i<=MAXTIMES; i++)
   totalBytes += writeFile(datum, count);
  tIOEnd = MPI_Wtime();
  tION = (tIOEnd - tIOStart)/MAXTIMES;
  MPI_File_close (&fileHandle);

	/*
	 * * * * * * * * * Independent MPI-IO to file system from all compute nodes - shared file * * * * * * * *
	 */
  MPI_File_open (MPI_COMM_WORLD, fileNameFS, mode, MPI_INFO_NULL, &fileHandle);
	for (int i=1; i<=SKIP; i++)
	 totalBytes += writeFile(datum, count);
	tIOStart = MPI_Wtime();
	for (int i=1; i<=MAXTIMESD; i++)
	 totalBytes += writeFile(datum, count);
	tIOEnd = MPI_Wtime();
  tFS = (tIOEnd - tIOStart)/MAXTIMESD;
	MPI_File_close (&fileHandle);

#else

  MPI_File_open (MPI_COMM_WORLD, testFileName, mode, MPI_INFO_NULL, &fileHandle);
	for (int i=1; i<=SKIP; i++)
	 totalBytes += writeFile(datum, count);
	tIOStart = MPI_Wtime();
	for (int i=1; i<=MAXTIMESD; i++)
	 totalBytes += writeFile(datum, count);
	tIOEnd = MPI_Wtime();
  tFS = (tIOEnd - tIOStart)/MAXTIMESD;
	MPI_File_close (&fileHandle);

#endif

}

int readFile(dataBlock *datum, int count) 
{
  int result, nbytes;
  int blocking = 0;

	MPI_Request request;

	if (blocking == 1) {
		result = MPI_File_read_at (fileHandle, (MPI_Offset)rank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &status);
		if (result != MPI_SUCCESS) 
			prnerror (result, "all MPI_File_write_at Error:");
	}
	else {
		MPIO_Request req_iw;
		MPI_Status st_iw;
		result = MPI_File_iread_at (fileHandle, (MPI_Offset)rank*count*sizeof(double), datum->getAlphaBuffer(), count, MPI_DOUBLE, &req_iw);
		MPI_Wait(&req_iw, &st_iw);
	}

	MPI_Get_elements( &status, MPI_CHAR, &nbytes );

	return nbytes;
}

void file_read(dataBlock *datum) {

  int totalBytes = 0;
  double tIOStart, tIOEnd;

#ifdef BGQ
   // something
  puts(" ");
#else
  MPI_File_open (MPI_COMM_WORLD, testFileName, mode, MPI_INFO_NULL, &fileHandle);
	for (int i=1; i<=SKIP; i++)
	 totalBytes += readFile(datum, count);
	tIOStart = MPI_Wtime();
	for (int i=1; i<=MAXTIMESD; i++)
	 totalBytes += readFile(datum, count);
	tIOEnd = MPI_Wtime();
  tFS = (tIOEnd - tIOStart)/MAXTIMESD;
	MPI_File_close (&fileHandle);
#endif  

}

int main(int argc, char *argv[]) {

  int max_tION, max_tFS;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  sprintf(testFileName,"TestFile-%d",size);

#ifndef THETA
  init_(argc, argv);
#endif
  count = atoi(argv[1]) * oneKB;

	/* allocate buffer */
  dataBlock *datum = new dataBlock(count);
	datum->allocElement ();

  assert (datum->getAlphaBuffer() != NULL);

  /* set file open mode */
  mode = MPI_MODE_CREATE | MPI_MODE_RDWR; //WRONLY;

  file_write(datum);
  file_read(datum);
  
	MPI_Reduce(&tION, &max_tION, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tFS, &max_tFS, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (rank == 0)
  printf ("0: Times: %d | %6.2f KB | %4.2lf %4.2lf\n", size, 8.0*count/1024.0, max_tION, max_tFS);

  MPI_Finalize();

}

