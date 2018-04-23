#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#define RUNS 100

// Macros used in reduce collective

int main(int argc,char *argv[]){
	int rank, p;
	char* ptr;
    MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 
	if(argc!=2){	 
		if(rank==0) printf("Usage: <program> <message_size>\n");
		exit(0);
	}	
	long int SIZE = strtol(argv[1], &ptr, 10);
	int len = SIZE / sizeof(int); 

	int *recvmsg = malloc(len * sizeof(int));
	int *selfmsg = malloc(len * sizeof(int));
	int *reducedmsg = malloc(len * sizeof(int));
	
	double t1,t2,res;

  // Do the recusrive halfing
	for (int i=0;i<RUNS;i++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		for (int ll = 0; ll < len; ll++) {
			selfmsg[ll] = ll;
			reducedmsg[ll] = ll;
		}	

		t1 = MPI_Wtime();
		  int vrank = rank;
			int xy = floor(log2(p));
		  unsigned pof2 = pow(2, xy);
		  int rem = p - pof2;
		  unsigned mask = 1;

		  if (rank < 2 * rem) {
			if(rank % 2 == 0) {

			  MPI_Send(selfmsg, len, MPI_INT, rank+1, 1, MPI_COMM_WORLD);
			  vrank = -1;
			} else {
			  
			  MPI_Recv(recvmsg, len, MPI_INT, rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			 
				for (int ll = 0; ll < len; ll++) {
					reducedmsg[ll] += recvmsg[ll];
				}
			  vrank = rank / 2;
			}
		  } else {
			vrank = rank - rem;
		  }

		  if(vrank != -1) {
			// Recursive doubling
			while(mask < pof2) {
			  int dest = vrank ^ mask; // bitwise xor
			  if(dest < rem) {
				dest = dest*2 + 1;
			  } else {
				dest = dest + rem;
			  }

			  MPI_Sendrecv(reducedmsg, len, MPI_INT, dest, 1, recvmsg, len, MPI_INT, dest, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			  
			  for (int ll = 0; ll < len; ll++) {
					reducedmsg[ll] += recvmsg[ll];
				}
			  
			  mask = mask << 1;
			}
		  }

		  if (rank < 2*rem) {
			// Expand from power of two
			if(rank % 2 > 0) {
			  MPI_Send(reducedmsg, len, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
			} else {
			  MPI_Recv(reducedmsg, len, MPI_INT, rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			  
			}
		  }
			
		t2 = MPI_Wtime() - t1;
		MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		if(rank==0)
		{
			printf("Run %ld time %1.9lf\n", i+1,res);

		}
	}
 	MPI_Finalize();

}
