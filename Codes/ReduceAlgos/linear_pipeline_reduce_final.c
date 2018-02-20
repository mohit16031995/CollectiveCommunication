#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define RUNS 10

//#define SIZE 6
//#define CHUNK 2
//#define CSIZE SIZE/CHUNK

int main(int argc, char *argv[]) {

	int rank, p, root = 0, index, cdone = 0, leaf, k ,current_chunk_start = 0;
	int peer;
	long int count, i, j, SIZE, CSIZE;
	char *ptr;
	int CHUNK;
	int recvFromChild = -1;
	int children = 0;
	int curCh = 0;
	int recv = 1;
	
	int seed = time(NULL);
	srand(seed);

	//char inmsg[SIZE], selfmsg[SIZE];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	leaf = p - 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc != 3) {
		if (rank == 0) 
			printf("Usage: <program> <message_size> <nchunk>\n");
		exit(0);
	}

	SIZE = strtol (argv[1], &ptr, 10); //length of int vector
	CHUNK = atoi (argv[2]); //no of chunks
	CSIZE = SIZE / CHUNK; //no of elements in each chunk
	SIZE = CSIZE * CHUNK; //total size augmented


	int *msg = malloc(SIZE * sizeof(int));
	int *selfmsg = malloc(SIZE * sizeof(int));

	// for recvs
	MPI_Status stt;
	MPI_Request req[CHUNK];

	// for send
	MPI_Status sstt[CHUNK*2];
	MPI_Request sreq[CHUNK*2];

	double t1, t2, res;

	for (i = 0; i < SIZE; i++) {
		selfmsg[i] = (i + rank) % 100;
		/*if(rank==root)*/ 
		msg[i] = selfmsg[i];
	}

	//selfmsg[SIZE] = 100;
	//msg[SIZE] = 100;

	if (rank != leaf) {
		recvFromChild = rank + 1;
	}

	if (rank != root) {
		children = 1;
                peer = rank - 1;		
	}

	double timings[2][50][515];

	for (i = 0; i < RUNS; i++) {

		MPI_Barrier(MPI_COMM_WORLD);
		cdone = 0;
		count = 0;
		t1 = MPI_Wtime();

		// set up all recv from left and right tree
		if (rank != leaf) {
			// if not leaf setup all recvs
			for (j = 0; j < CHUNK; j++) {
				MPI_Irecv(msg + j * CSIZE, CSIZE, MPI_INT, recvFromChild, j, MPI_COMM_WORLD, &req[j]);
			}
		} else {  // if leaf then setup all send's
			
			for(j = 0; j < CHUNK; j++) {
				
				if(children) {
					MPI_Isend(selfmsg + j * CSIZE, CSIZE, MPI_INT, peer, j, MPI_COMM_WORLD, &sreq[count++]);
					timings[1][j][peer] = MPI_Wtime() - t1;
					
				}
			}
		}

		//if not leaf -> check which recv finish and setup send for them
		if(rank != leaf)
			while(cdone < CHUNK) {
				MPI_Waitany(CHUNK, req, &index, &stt);
				if(index == MPI_UNDEFINED) {
					printf("Unexpected error!\n");
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				current_chunk_start = index*CSIZE;
				for (k = current_chunk_start; k < current_chunk_start + CSIZE; k++)
				{
					msg[k] = (msg[k] + selfmsg[k])%100;
				}

				timings[0][index][recvFromChild] = MPI_Wtime() - t1;
				if(children) {
						MPI_Isend(msg + index * CSIZE, CSIZE, MPI_INT, peer, index, MPI_COMM_WORLD, &sreq[count++]);
						timings[1][index][peer] = MPI_Wtime() - t1;
				}
				cdone++;
			}

		MPI_Waitall(count, sreq, sstt);  // wait for all send to finish
		t2 = MPI_Wtime() - t1;
		MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		if(rank == 0){
			printf("Run %ld time %1.9lf\n\n", i + 1, res);
			for (int ii = 0; ii < SIZE; ii++) {
				printf ("%d  ",msg[ii]);
			}
			printf("\n");
		} else {;
			//printf("selfmsg %s Inmsg %s\n", selfmsg,msg);
			//j=strcmp(selfmsg,msg);
			//j!=0?printf("Error - msgs different\n"):0;
			//memset(msg,'$',SIZE);
		}		
	}
	MPI_Finalize();
}
