#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//#define RUNS 100
//#define SIZE 6
//#define CHUNK 2
//#define CSIZE SIZE/CHUNK

// Macros used in reduce collective

int main(int argc,char *argv[]){
	int rank,p,root = 0, index,cdone=0;
	int vrank, vpeer,peer;
	long int i, count, SIZE,CSIZE;
	char *ptr;

	int CHUNK;

	int recvToLeft = -1;
    	int recvToRight = -1;

	int leftChildren = 0;
	int rightChildren = 0;

 	int leftPeers[2];
	int rightPeers[2];


//	int seed = time(NULL);
//	srand(seed);
	//char inmsg[SIZE], outmsg[SIZE];
//
//	MPI_Init(&argc,&argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &p);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(argc!=4){	 
	  if(rank==0) printf("Usage: <program> <message_size> <nchunk>\n");
	  exit(0);
	}	
	
	int nP = strtol(argv[1], &ptr, 10);
	SIZE = strtol(argv[2], &ptr, 10);
	CHUNK = atoi(argv[3]);
	CSIZE = SIZE/CHUNK; 
	SIZE = CSIZE*CHUNK;

//	char *msg = malloc(SIZE+1 * sizeof(char));
//	char *outmsg = malloc(SIZE+1 * sizeof(char));

	// for recvs
//	MPI_Status stt;
//	MPI_Request req[CHUNK];

	// for send
//	MPI_Status sstt[CHUNK*2];
//	MPI_Request sreq[CHUNK*2];
	
//	double t1,t2,res;

//	for(i=0;i<SIZE;i++) {
//		outmsg[i] = 'A'+i%26;
//		if(rank==0) msg[i] = outmsg[i];	
//	}
//	outmsg[SIZE] = '\0';	
//	msg[SIZE] = '\0';
	
	printf("num_ranks %d\n", nP);
	p = nP;
	count = 1;
	printf("\nrank 0 {\n");
	long int time = 0;

	///////////////////////////////////////////////////milliseconds
	time = ((3*CSIZE*100)/(sizeof(int)));
	///////////////////////////////////////////////////

	for (i = 0; i < CHUNK; i++) {
		if (i % 2 == 0) {
				printf("l%ld: recv %ldb from 1 tag %ld\n", count++, CSIZE, i+1);
				//printf("l%ld: send %ldb to 1\n", count++, CSIZE);
		}
		else {
				printf("l%ld: recv %ldb from %d tag %ld\n", count++, CSIZE, p-1, i+1);
				//printf("l%ld: send %ldb to %d\n", count++, CSIZE, p-1);
		}
	}	
	for (i = 0; i < CHUNK; i++) {
		printf("l%ld: calc %ld\n", count, time);
		printf("l%ld requires l%ld\n", count++, i+1);
	}
	printf("}\n");
	for (rank = 1; rank < nP; rank++){
		count = 1;
		leftChildren = 0;
		rightChildren = 0;
		recvToLeft = rank / 2;
		recvToRight = (p - (int)floor((p-rank)/2))%p;
		printf("\nrank %d {\n", rank);
		/*for (i = 0; i < CHUNK; i++) {
			if (i % 2 == 0) {
				printf("l%ld: recv %ldb from %d tag %ld\n", count++, CSIZE, recvToLeft, i+1);
				//printf("l%ld: recv %ldb from %d\n", count++, CSIZE, recvToLeft);
			}
			else {
				printf("l%ld: recv %ldb from %d tag %ld\n", count++, CSIZE, recvToRight, i+1);
				//printf("l%ld: recv %ldb from %d\n", count++, CSIZE, recvToRight);
			}
		}*/	
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
		for (i = 0; i < CHUNK; i++)  {
			if (i % 2 == 0 && leftChildren != 0) {
				printf("l%ld: recv %ldb from %d tag %ld\n", count++, CSIZE, leftPeers[0], i+1);
				//printf("l%ld: send %ldb to %d\nl%ld irequires l%ld\n", count, CSIZE, leftPeers[0], count, i+1);
				if (leftChildren == 2) {
					printf("l%ld: recv %ldb from %d tag %ld\n", count++, CSIZE, leftPeers[1], i+1);
					//printf("l%ld: send %ldb to %d\nl%ld irequires l%ld\n", count, CSIZE, leftPeers[1], count, i+1);	
				}
			}
			else if (i % 2 != 0 && rightChildren != 0) {
				printf("l%ld: recv %ldb from %d tag %ld\n", count++, CSIZE, rightPeers[0], i+1);
				//printf("l%ld: send %ldb to %d\nl%ld irequires l%ld\n", count, CSIZE, rightPeers[0], count, i+1);
				if (rightChildren == 2) {
					printf("l%ld: recv %ldb from %d tag %ld\n", count++, CSIZE, rightPeers[1], i+1);
					//printf("l%ld: send %ldb to %d\nl%ld irequires l%ld\n", count, CSIZE, rightPeers[0], count, i+1);
				}
			}
		}
		for (i = 0; i < CHUNK; i++) {
			if (i % 2 == 0) {
				if (leftChildren == 0) {
					printf("l%ld: send %ldb to %d tag %ld\n", count++, CSIZE, recvToLeft, i+1);
				}
			}
			else {
				if (rightChildren == 0) {
					printf("l%ld: send %ldb to %d tag %ld\n", count++, CSIZE, recvToRight, i+1);
				}
			}
		}
		int count2 = 1;
		for (i = 0; i < CHUNK; i++) {
			if (i % 2 == 0) {
				if (leftChildren == 1) {
					printf("l%d: calc %lld\n", count, time);
					printf("l%d requires l%d\n", count, count2++);
					printf("l%ld: send %ldb to %d tag %ld\n", count+1, CSIZE, recvToLeft, i+1);
					printf("l%d requires l%d\n", count+1, count);
					count += 2;
				}
				else if (leftChildren == 2) {
					printf("l%d: calc %lld\n", count, time);
					printf("l%d requires l%d\n", count, count2++);
					printf("l%d: calc %lld\n", count+1, time);
					printf("l%d requires l%d\n", count+1, count2++);
					printf("l%ld: send %ldb to %d tag %ld\n", count+2, CSIZE, recvToLeft, i+1);
					printf("l%d requires l%d\n", count+2, count);
					printf("l%d requires l%d\n", count+2, count+1);
					count += 3;
				}
			}
			else {
				if (rightChildren == 1) {
					printf("l%d: calc %lld\n", count, time);
					printf("l%d requires l%d\n", count, count2++);
					printf("l%ld: send %ldb to %d tag %ld\n", count+1, CSIZE, recvToRight, i+1);
					printf("l%d requires l%d\n", count+1, count);
					count += 2;
				}
				else if (rightChildren == 2){
					printf("l%d: calc %lld\n", count, time);
					printf("l%d requires l%d\n", count, count2++);
					printf("l%d: calc %lld\n", count+1, time);
					printf("l%d requires l%d\n", count+1, count2++);
					printf("l%ld: send %ldb to %d tag %ld\n", count+2, CSIZE, recvToRight, i+1);
					printf("l%d requires l%d\n", count+2, count);
					printf("l%d requires l%d\n", count+2, count+1);
					count += 3;
				}
			}
		}	
		printf("}\n");
	}
	return 0;
}
