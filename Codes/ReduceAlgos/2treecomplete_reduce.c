#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define RUNS 100
//#define SIZE 6
//#define CHUNK 2
//#define CSIZE SIZE/CHUNK

// Macros used in reduce collective

int main(int argc,char *argv[]){
	int rank, p, root = 0, index, cdone=0;
	long int count, i, j, k, SIZE, CSIZE, pseudo_index;
	int CHUNK;

	int parentLeft = -1;
    	int parentRight = -1;

	int leftChildren = 1;
	int rightChildren = 1;

 	int leftPeers[2] = {1,-1};
	int rightPeers[2];


	int seed = time(NULL);
	srand(seed);
	//char inmsg[SIZE], selfmsg[SIZE];

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	rightPeers[0] = p-1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(argc!=3){	 
	  if(rank==0) printf("Usage: <program> <message_size> <nchunk>\n");
	  exit(0);
	}	
	
	SIZE = strtol(argv[1], &ptr, 10);
	CHUNK = atoi(argv[2]);
	CSIZE = SIZE/CHUNK; 
	SIZE = CSIZE*CHUNK;

	int *msg1 = malloc(SIZE*2 * sizeof(int));
	int *selfmsg = malloc(SIZE * sizeof(int));
	int *msg2 = msg1+CHUNK;
	int ready[CHUNK] = {0};

	// for  recvs 
	MPI_Status stt;
	MPI_Request req1[CHUNK*2];
	MPI_Request *req2 = req1+CHUNK;

	// for  send
	MPI_Status sstt[CHUNK];
	MPI_Request sreq[CHUNK];
	
	double t1,t2,res;

	for (i=0;i<SIZE;i++) {
		selfmsg[i] = (i%100 + p%100)%100;
		//if(rank==0) msg[i] = selfmsg[i];	
	}
	//selfmsg[SIZE] = '\0';	
	//msg[SIZE] = '\0';
	
	
	if (rank != 0) 						//if not root
	{
		parentLeft = rank / 2;
		parentRight = (p - (p-rank)/2)%p;
		if (2*rank < p) {
			leftChildren = 1;
			leftPeers[0] = (2*rank);
		}
		if ((2*rank)+1 < p) {
			leftChildren = 2;
			leftPeers[1] = (2*rank)+1;
		}
		if (2*rank - p > 0) {
			rightChildren = 1;
			rightPeers[0] = (2*rank - p);
		}
		if (2*rank-p-1 > 0) {
			rightChildren = 2;
			rightPeers[1] = (2*rank-p-1);
		}
	}

//	double timings[2][50][515];

	for (i=0;i<RUNS;i++)
	{
		MPI_Barrier(MPI_COMM_WORLD);

		cdone=0; count=0;
		for (k=0;k<CHUNK;k++)
			ready[k]=0;
		t1 = MPI_Wtime();
		
	// set up all recv from left and right tree
		if (leftChildren)				//if not leafleft setup all even recvs
			for (j=0;j<CHUNK;j+=2) {
				MPI_Irecv(msg1+j*CSIZE,CSIZE,MPI_INT,leftPeers[0],j,MPI_COMM_WORLD,&req1[j]);
				if (leftChildren==2)
				MPI_Irecv(msg2+j*CSIZE,CSIZE,MPI_INT,leftPeers[1],j,MPI_COMM_WORLD,&req2[j]);
			}
		if (rightChildren)				//if not leafRight setup all odd recvs
			for (j=1;j<CHUNK;j+=2) {
				MPI_Irecv(msg1+j*CSIZE,CSIZE,MPI_INT,rightPeers[0],j,MPI_COMM_WORLD,&req1[j]);
				if (rightChildren==2)
				MPI_Irecv(msg2+j*CSIZE,CSIZE,MPI_INT,rightPeers[1],j,MPI_COMM_WORLD,&req2[j]);
			}

	// setup all sends in left and right leaves
		if (!leftChildren || !rightChildren)
		for (j=0;j<CHUNK;j++) {	
			if (!leftChildren && !(j%2)) 
				MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,parentLeft,j,MPI_COMM_WORLD,&sreq[j]);
			if (!rightChildren && (j%2))
				MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,parentRight,j,MPI_COMM_WORLD,&sreq[j]);
		}

	// check which recv finish and setup send for them
		int chunks_to_recv = ((CHUNKS+1)/2)*leftChildren + (CHUNK/2)*rightChildren;
		if(rank!=root) 
			while(cdone < (chunks_to_recv) {		
				MPI_Waitany(CHUNK*2, req1, &index, &stt);
				if(index == MPI_UNDEFINED) {
					printf("Unexpected error!\n");
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				pseudo_index = (index>=CHUNK) ? (index-CHUNK)*CSIZE : index*CSIZE;
				for (k=pseudo_index,j=index*CSIZE;k<pseudo_index+CSIZE;k++,j++) {
					selfmsg[k] = (msg1[j]%100+selfmsg[k]%100)%100;
				}
				if ( (pseudo_index%2) ? ++ready[pseudo_index] == rightChildren : ++ready[pseudo_index] == leftChildren )
				{	//forward this chunk ;
					if (!pseudo_index%2)
						MPI_Isend(selfmsg+pseudo_index,CSIZE,MPI_INT,parentLeft,pseudo_index,MPI_COMM_WORLD,&sreq[count++]);
					else
						MPI_Isend(selfmsg+pseudo_index,CSIZE,MPI_INT,parentRight,pseudo_index,MPI_COMM_WORLD,&sreq[count++]);
				}
				cdone++;
			  }
		else
			while(cdone < (chunks_to_recv) {		
				MPI_Waitany(CHUNK*2, req1, &index, &stt);
				if(index == MPI_UNDEFINED) {
					printf("Unexpected error!\n");
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				pseudo_index = (index>=CHUNK) ? (index-CHUNK)*CSIZE : index*CSIZE;
				for (k=pseudo_index,j=index*CSIZE;k<pseudo_index+CSIZE;k++,j++) {
					selfmsg[k] = (msg1[j]%100+selfmsg[k]%100)%100;
				}
				cdone++;
			  }
		


		MPI_Waitall(count,sreq,sstt);  // wait for  all send to finish
		
		t2 = MPI_Wtime() - t1;
		MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		if(rank==0){
			printf("Run %ld time %1.9lf\n", i+1,res);
		} else {
			//printf("selfmsg %s Inmsg %s\n", selfmsg,msg);
			j=strcmp(selfmsg,msg);
			j!=0?printf("Error - msgs different\n"):0;
			memset(msg,'$',SIZE);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if (rank != 0) {
			for  (int ll = 0; ll < CHUNK; ll++) {
				int from;
				if (ll % 2) {
					from = parentRight;
					if (rightChildren==2) {
						printf("Logs, Process %d, Run %ld, chunk %d, received %d %1.9lf, left_sent %d %1.9lf, right_sent %d %1.9lf\n", rank, i+1, ll, from, timings[0][ll][from], rightPeers[0], timings[1][ll][rightPeers[0]], rightPeers[1], timings[1][ll][rightPeers[1]]);
					}
					else if (rightChildren) {
						printf("Logs, Process %d, Run %ld, chunk %d, received %d %1.9lf, left_sent %d %1.9lf\n", rank, i+1, ll, from, timings[0][ll][from], rightPeers[0], timings[1][ll][rightPeers[0]]);
					}
					else {
						printf("Logs, Process %d, Run %ld, chunk %d, received %d %1.9lf\n", rank, i+1, ll, from, timings[0][ll][from]);
					}
				}
				else {
					from = parentLeft;
					if (leftChildren == 2) {
						printf("Logs, Process %d, Run %ld, chunk %d, received %d %1.9lf, left_sent %d %1.9lf, right_sent %d %1.9lf\n", rank, i+1, ll, from, timings[0][ll][from], leftPeers[0], timings[1][ll][leftPeers[0]], leftPeers[1], timings[1][ll][leftPeers[1]]);
					}	
					else if (leftChildren) {
						printf("Logs, Process %d, Run %ld, chunk %d, received %d %1.9lf, left_sent %d %1.9lf\n", rank, i+1, ll, from, timings[0][ll][from], leftPeers[0], timings[1][ll][leftPeers[0]]);
					}	
					else {
						printf("Logs, Process %d, Run %ld, chunk %d, received %d %1.9lf\n", rank, i+1, ll, from, timings[0][ll][from]);
					}			
				}
			}
		}
		else {
			for (int ll=0;ll<CHUNK;ll++) {
			// left tree				
			if(!(ll%2)) { 
				printf("Logs, Process %d, Run %ld, chunk %d, left_sent %d %1.9lf\n", rank, i+1, ll, 1, timings[1][ll][1]);
			}			
			// right tree
			else { 
				printf("Logs, Process %d, Run %ld, chunk %d, left_sent %d %1.9lf\n", rank, i+1, ll, p-1, timings[1][ll][p-1]);
			}			
		  	}
		}
	}
	MPI_Finalize();
} 

