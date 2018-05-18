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


int main(int argc,char *argv[]){
	int rank,p,root=0,index,cdone=0,leaf,k;
	int vrank, vpeer,peer;
	long int count,i,j,SIZE,CSIZE;
	char *ptr;

	int CHUNK;

	int recvToLeft = -1;
    	int recvToRight = -1;

	int leftChildren = 0;
	int rightChildren = 0;

 	int leftPeers[2];
	int rightPeers[2];

	int curChLeft = 0;
	int curChRight = 1;

	int recvLeft = 1;
	int recvRight = 1;


	int seed = time(NULL);
	srand(seed);
	//char inmsg[SIZE], selfmsg[SIZE];

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	leaf = p-1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(argc!=3){
	  if(rank==0) printf("Usage: <program> <message_size> <nchunk>\n");
	  exit(0);
	}

	SIZE = strtol(argv[1], &ptr, 10); //length of int vector
	CHUNK = atoi(argv[2]); //no of chunks
	CSIZE = SIZE/CHUNK; //no of elements in each chunk
	SIZE = CSIZE*CHUNK; //total size augmented


	int *msg = malloc(SIZE * sizeof(int));
	int *selfmsg = malloc(SIZE * sizeof(int));

	// for recvs
	MPI_Status stt;
	MPI_Request req[CHUNK];

	// for send
	MPI_Status sstt[CHUNK*2];
	MPI_Request sreq[CHUNK*2];

	double t1,t2,res;

	for(i=0;i<SIZE;i++) {
		selfmsg[i] = (i+rank)%100;
		/*if(rank==root)*/ msg[i] = selfmsg[i];
	}
	//selfmsg[SIZE] = 100;
	//msg[SIZE] = 100;

	if (rank != leaf) {
		recvToLeft = rank+1;
		recvToRight = rank+1;
	}
	if (rank != root) {
		leftChildren = 1;
		rightChildren = 1;
		leftPeers[0] = rank-1;
		rightPeers[0] = rank-1;
	}


	double timings[2][50][515];
	for(i=0;i<RUNS;i++){
		MPI_Barrier(MPI_COMM_WORLD);

		cdone=0; count=0;
		t1 = MPI_Wtime();

		// set up all recv from left and right tree
		if(rank!=leaf) {// if not leaf setup all recvs
		   for(j=0;j<CHUNK;j++) {
			// left tree
			if(!(j%2)) MPI_Irecv(msg+j*CSIZE,CSIZE,MPI_INT,recvToLeft,j,MPI_COMM_WORLD,&req[j]);
			// right tree
			else MPI_Irecv(msg+j*CSIZE,CSIZE,MPI_INT,recvToRight,j,MPI_COMM_WORLD,&req[j]);
		   }
		} else {  // if leaf then setup all send's
		   for(j=0;j<CHUNK;j++) {
			// left tree
			if(!(j%2)) {
				if(leftChildren) {
					MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,leftPeers[0],j,MPI_COMM_WORLD,&sreq[count++]);
					timings[1][j][leftPeers[0]] = MPI_Wtime() - t1;
				}
				if(leftChildren==2) {
					MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,leftPeers[1],j,MPI_COMM_WORLD,&sreq[count++]);
					timings[1][j][leftPeers[1]] = MPI_Wtime() - t1;
				}
			}
			// right tree
			else {
				if(rightChildren) {
					MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,rightPeers[0],j,MPI_COMM_WORLD,&sreq[count++]);
					timings[1][j][rightPeers[0]] = MPI_Wtime() - t1;
				}
				if(rightChildren==2) {
					MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,rightPeers[1],j,MPI_COMM_WORLD,&sreq[count++]);
					timings[1][j][rightPeers[1]] = MPI_Wtime() - t1;
				}
			}
		   }
		}

		//if not leaf -> check which recv finish and setup send for them
		if(rank!=leaf)
		  while(cdone < CHUNK) {
			MPI_Waitany(CHUNK, req, &index, &stt);
	            	if(index == MPI_UNDEFINED) {
                		printf("Unexpected error!\n");
                		MPI_Abort(MPI_COMM_WORLD, 1);
            		}

								for (k=current_chunk_start; k<current_chunk_start+CSIZE; k++)
								{
										msg[k] = (msg[k] + selfmsg[k])%100;
								}


			if(index%2) {  // right tree
				timings[0][index][recvToRight] = MPI_Wtime() - t1;
				if(rightChildren) {
					current_chunk_start = msg+index*CSIZE;
					MPI_Isend(msg+index*CSIZE,CSIZE,MPI_INT,rightPeers[0],index,MPI_COMM_WORLD,&sreq[count++]);
					timings[1][index][rightPeers[0]] = MPI_Wtime() - t1;
				}
				if(rightChildren==2) {
					MPI_Isend(msg+index*CSIZE,CSIZE,MPI_INT,rightPeers[1],index,MPI_COMM_WORLD,&sreq[count++]);
					timings[1][index][rightPeers[1]] = MPI_Wtime() - t1;
				}
			} else {
				timings[0][index][recvToLeft] = MPI_Wtime() - t1;
				if(leftChildren) {
					MPI_Isend(msg+index*CSIZE,CSIZE,MPI_INT,leftPeers[0],index,MPI_COMM_WORLD,&sreq[count++]);
					timings[1][index][leftPeers[0]] = MPI_Wtime() - t1;
				}
				if(leftChildren==2) {
					MPI_Isend(msg+index*CSIZE,CSIZE,MPI_INT,leftPeers[1],index,MPI_COMM_WORLD,&sreq[count++]);
					timings[1][index][leftPeers[1]] = MPI_Wtime() - t1;
				}
			}
			cdone++;
		  }


		MPI_Waitall(count,sreq,sstt);  // wait for all send to finish

		t2 = MPI_Wtime() - t1;
		MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		if(rank==0){
			printf("Run %ld time %1.9lf\n\n", i+1,res);

			for (i=0; i<SIZE; i++)
			{
					printf ("%d  ",msg[i]);
			}
			printf("\n");
		} else {;
			//printf("selfmsg %s Inmsg %s\n", selfmsg,msg);
			//j=strcmp(selfmsg,msg);
			//j!=0?printf("Error - msgs different\n"):0;
			//memset(msg,'$',SIZE);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if (rank != 0) {
			for (int ll = 0; ll < CHUNK; ll++) {
				int from;
				if (ll % 2) {
					from = recvToRight;
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
					from = recvToLeft;
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
			for(int ll=0;ll<CHUNK;ll++) {
			// left tree
			if(!(ll%2)) {
				if(leftChildren==2) {
					printf("Logs, Process %d, Run %ld, chunk %d, left_sent %d %1.9lf, right_sent %d %1.9lf\n", rank, i+1, ll, leftPeers[0], timings[1][ll][leftPeers[0]], leftPeers[1], timings[1][ll][leftPeers[1]]);
				}
				else if(leftChildren) {
					printf("Logs, Process %d, Run %ld, chunk %d, left_sent %d %1.9lf\n", rank, i+1, ll, leftPeers[0], timings[1][ll][leftPeers[0]]);
				}
			}
			// right tree
			else {
				if(rightChildren==2) {
					printf("Logs, Process %d, Run %ld, chunk %d, left_sent %d %1.9lf, right_sent %d %1.9lf\n", rank, i+1, ll, rightPeers[0], timings[1][ll][rightPeers[0]], rightPeers[1], timings[1][ll][rightPeers[1]]);
				}
				else if(rightChildren) {
					printf("Logs, Process %d, Run %ld, chunk %d, left_sent %d %1.9lf\n", rank, i+1, ll, rightPeers[0], timings[1][ll][rightPeers[0]]);
				}
			}
		   	}
		}
	}
	MPI_Finalize();
}
