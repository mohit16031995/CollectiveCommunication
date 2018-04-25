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
	for (i = 0; i < CHUNK; i++) {
	
				printf("l%ld: send %ldb to 1 tag %ld\n", count++, CSIZE, i+1);
				//printf("l%ld: send %ldb to 1\n", count++, CSIZE);
		
	
				printf("l%ld: send %ldb to 2 tag %ld\n", count++, CSIZE,i+1);
				//printf("l%ld: send %ldb to %d\n", count++, CSIZE, p-1);

	}	
	printf("}\n");
	for (rank = 1; rank < nP; rank++){
		count = 1;
		leftChildren = 0;
		rightChildren = 0;
		recvToLeft = ((rank+1)/2)-1;
		recvToRight = ((rank+1)/2)-1;
		printf("\nrank %d {\n", rank);
		for (i = 0; i < CHUNK; i++) {
			if (i % 2 == 0) {
				printf("l%ld: recv %ldb from %d tag %ld\n", count++, CSIZE, recvToLeft, i+1);
				//printf("l%ld: recv %ldb from %d\n", count++, CSIZE, recvToLeft);
			}
			else {
				printf("l%ld: recv %ldb from %d tag %ld\n", count++, CSIZE, recvToRight, i+1);
				//printf("l%ld: recv %ldb from %d\n", count++, CSIZE, recvToRight);
			}
		}	
		if ((2*rank)+1 < p) {
			leftChildren = 1;
			leftPeers[0] = (2*rank)+1;
			rightChildren = 1;
			rightPeers[0] = (2*rank)+1;
		}
		if ((2*rank)+2 < p) {
			leftChildren = 2;
			leftPeers[1] = (2*rank)+2;
			rightChildren = 2;
			rightPeers[1] = (2*rank)+2;
		}

		for (i = 0; i < CHUNK; i++)  {
			if (i % 2 == 0 && leftChildren != 0) {
				printf("l%ld: send %ldb to %d tag %ld\nl%ld requires l%ld\n", count, CSIZE, leftPeers[0], i+1, count, i+1);
				//printf("l%ld: send %ldb to %d\nl%ld irequires l%ld\n", count, CSIZE, leftPeers[0], count, i+1);
				count++;
				if (leftChildren == 2) {
					printf("l%ld: send %ldb to %d tag %ld\nl%ld requires l%ld\n", count, CSIZE, leftPeers[1], i+1, count, i+1);
					//printf("l%ld: send %ldb to %d\nl%ld irequires l%ld\n", count, CSIZE, leftPeers[1], count, i+1);
					count++;	
				}
			}
			else if (i % 2 != 0 && rightChildren != 0) {
				printf("l%ld: send %ldb to %d tag %ld\nl%ld requires l%ld\n", count, CSIZE, rightPeers[0], i+1, count, i+1);
				//printf("l%ld: send %ldb to %d\nl%ld irequires l%ld\n", count, CSIZE, rightPeers[0], count, i+1);
				count++;
				if (rightChildren == 2) {
					printf("l%ld: send %ldb to %d tag %ld\nl%ld requires l%ld\n", count, CSIZE, rightPeers[1], i+1, count, i+1);
					//printf("l%ld: send %ldb to %d\nl%ld irequires l%ld\n", count, CSIZE, rightPeers[0], count, i+1);
					count++;
				}
			}
		}
		printf("}\n");
	}
	return 0;
}
