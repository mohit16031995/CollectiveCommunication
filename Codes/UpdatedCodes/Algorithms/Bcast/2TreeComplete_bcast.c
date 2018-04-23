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
	int rank,p,root = 0, index,cdone=0;
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


	int seed = time(NULL);
	srand(seed);
	//char inmsg[SIZE], outmsg[SIZE];

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(argc!=3){	 
	  if(rank==0) printf("Usage: <program> <message_size> <nchunk>\n");
	  exit(0);
	}	
	
	SIZE = strtol(argv[1], &ptr, 10);
	CHUNK = atoi(argv[2]);
	CSIZE = SIZE/CHUNK; 
	SIZE = CSIZE*CHUNK;

	char *msg = malloc(SIZE+1 * sizeof(char));
	char *outmsg = malloc(SIZE+1 * sizeof(char));

	// for recvs
	MPI_Status stt;
	MPI_Request req[CHUNK];

	// for send
	MPI_Status sstt[CHUNK*2];
	MPI_Request sreq[CHUNK*2];
	
	double t1,t2,res;

	for(i=0;i<SIZE;i++) {
		outmsg[i] = 'A'+i%26;
		if(rank==0) msg[i] = outmsg[i];	
	}
	outmsg[SIZE] = '\0';	
	msg[SIZE] = '\0';
	
	
	if (rank != 0) {
		recvToLeft = rank / 2;
		recvToRight = (p - (int)floor((p-rank)/2))%p;
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
	for(i=0;i<RUNS;i++){
		MPI_Barrier(MPI_COMM_WORLD);

		cdone=0; count=0;
		t1 = MPI_Wtime();
		
		// set up all recv from left and right tree
		if(rank!=0) {// if not root setup all recvs
		   for(j=0;j<CHUNK;j++) {
			// left tree				
			if(!(j%2)) MPI_Irecv(msg+j*CSIZE,CSIZE,MPI_CHAR,recvToLeft,j,MPI_COMM_WORLD,&req[j]);			
			// right tree
			else MPI_Irecv(msg+j*CSIZE,CSIZE,MPI_CHAR,recvToRight,j,MPI_COMM_WORLD,&req[j]);
		   }
		} else {  // if root then setup all send's
		   for(j=0;j<CHUNK;j++) {
			// left tree				
			if(!(j%2)) { 
				MPI_Isend(outmsg+j*CSIZE,CSIZE,MPI_CHAR,1,j,MPI_COMM_WORLD,&sreq[count++]);
			}			
			// right tree
			else { 
				MPI_Isend(outmsg+j*CSIZE,CSIZE,MPI_CHAR,p-1,j,MPI_COMM_WORLD,&sreq[count++]);
			}			
		   }
		}						

		//if not root -> check which recv finish and setup send for them		
		if(rank!=root) 
		  while(cdone < CHUNK) {		
			MPI_Waitany(CHUNK, req, &index, &stt);
	            	if(index == MPI_UNDEFINED) {
                		printf("Unexpected error!\n");
                		MPI_Abort(MPI_COMM_WORLD, 1);
            		}
			if(index%2) {  // right tree
				if(rightChildren) {
					MPI_Isend(msg+index*CSIZE,CSIZE,MPI_CHAR,rightPeers[0],index,MPI_COMM_WORLD,&sreq[count++]);
				}
				if(rightChildren==2) {
					MPI_Isend(msg+index*CSIZE,CSIZE,MPI_CHAR,rightPeers[1],index,MPI_COMM_WORLD,&sreq[count++]);
				}
			} else {
				if(leftChildren) {
					MPI_Isend(msg+index*CSIZE,CSIZE,MPI_CHAR,leftPeers[0],index,MPI_COMM_WORLD,&sreq[count++]);
				}
				if(leftChildren==2) {
					MPI_Isend(msg+index*CSIZE,CSIZE,MPI_CHAR,leftPeers[1],index,MPI_COMM_WORLD,&sreq[count++]);
				}
			}
			cdone++;
		  }


		MPI_Waitall(count,sreq,sstt);  // wait for all send to finish
		
		t2 = MPI_Wtime() - t1;
		MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		if(rank==0){
			printf("Run %ld time %1.9lf\n", i+1,res);
		} else {
			//printf("Outmsg %s Inmsg %s\n", outmsg,msg);
			j=strcmp(outmsg,msg);
			j!=0?printf("Error - msgs different\n"):0;
			memset(msg,'$',SIZE);
		}
	}
	MPI_Finalize();
} 
