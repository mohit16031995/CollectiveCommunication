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
		
	int vrank = rank;
	int xy = floor(log2(p));
	unsigned pof2 = pow(2, xy);

	int len = SIZE / sizeof(int); 
	int clen = (len / pof2);
	int CSIZE = clen*sizeof(int);	
	SIZE = CSIZE*pof2;
	len = SIZE / sizeof(int);	


	int *recvmsg = malloc(len * sizeof(int));
	//int recvmsg[len];
	//int selfmsg[len];
	//int reducedmsg[len];
	int *selfmsg = malloc(len * sizeof(int));
	int *reducedmsg = malloc(len * sizeof(int));
	
	double t1,t2,res;
	//printf("rank");
	for (int i=0;i<RUNS;i++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		for (int ll = 0; ll < len; ll++) {
			selfmsg[ll] = ll+rank;
			reducedmsg[ll] = ll+rank;
		}	

		t1 = MPI_Wtime();
		int rem = p - pof2;
		

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
			// Recursive HALFING
			unsigned mask = pof2 >> 1;
			while(mask >= 1) {
			  int dest = vrank ^ mask; // bitwise xor
				int dest2;
			  if(dest < rem) {
				dest2 = dest*2 + 1;
			  } else {
				dest2 = dest + rem;
			  }
			  int send_start = (dest / mask);
				send_start *= mask;
			  int recv_start = (vrank / mask);
				recv_start *= mask;
			  int vlen = (len*mask) / pof2;
				//printf("rank %d send_start %d recv_start %d\n", rank, send_start, recv_start);
			  //MPI_Sendrecv(reducedmsg+(CSIZE*send_start), vlen, MPI_INT, dest2, 1, recvmsg+(CSIZE*recv_start), vlen, MPI_INT, dest2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&reducedmsg[clen*send_start], vlen, MPI_INT, dest2, 1, &recvmsg[clen*recv_start], vlen, MPI_INT, dest2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			  for (int ll = 0; ll < vlen; ll++) {
					reducedmsg[(clen*recv_start)+ll] += recvmsg[(clen*recv_start)+ll];
			  }
			  mask = mask >> 1;
			}
		}
		//printf("123\n");
		if(vrank != -1) {
			// Recursive doubling
			unsigned mask = 1;
			while(mask < pof2) {
			  int dest = vrank ^ mask; // bitwise xor
				int dest2;
			  if(dest < rem) {
				dest2 = dest*2 + 1;
			  } else {
				dest2 = dest + rem;
			  }
			  int send_start = (vrank / mask)*mask;
			  int recv_start = (dest / mask)*mask;
			 int vlen = (len*mask) / pof2;
			  //MPI_Sendrecv(reducedmsg+(CSIZE*send_start), clen, MPI_INT, dest2, 1, recvmsg+(CSIZE*recv_start), clen, MPI_INT, dest2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			  MPI_Sendrecv(&reducedmsg[clen*send_start], vlen, MPI_INT, dest2, 1, &recvmsg[clen*recv_start], vlen, MPI_INT, dest2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			  for (int ll = 0; ll < clen; ll++) {
					reducedmsg[(clen*recv_start)+ll] = recvmsg[(clen*recv_start)+ll];
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
