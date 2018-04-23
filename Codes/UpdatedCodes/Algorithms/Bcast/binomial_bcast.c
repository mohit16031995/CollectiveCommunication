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
	int rank,p,root=0,index,cdone=0;
	int vrank, vpeer,peer;
	long int count,i,j,SIZE,CSIZE;
	char *ptr;

	int CHUNK;

	
	int seed = time(NULL);
	srand(seed);
	//char inmsg[SIZE], outmsg[SIZE];

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(argc!=2){	 
	  if(rank==0) printf("Usage: <program> <message_size> <nchunk>\n");
	  exit(0);
	}	
	
	SIZE = strtol(argv[1], &ptr, 10);


	char *msg = malloc(SIZE+1 * sizeof(char));
	char *outmsg = malloc(SIZE+1 * sizeof(char));

	// for recvs
	MPI_Status stt;

	// for send
	MPI_Status sstt[500];
	MPI_Request sreq[500];
	
	double t1,t2,res;

	for(i=0;i<SIZE;i++) {
		outmsg[i] = 'A'+i%26;
		if(rank==root) msg[i] = outmsg[i];	
	}
	outmsg[SIZE] = '\0';	
	msg[SIZE] = '\0';
	int no_childs = 0;
	int childs[512];
	int parent;
	// ----- get send to and receive from on both trees -----
	 for(int r = 0; r < ceil(log2(p)); ++r) {
        int p2r = pow(2,r);

        int peer = rank + p2r;
        if(peer < p && rank < p2r) {
          childs[no_childs] = peer;
			no_childs += 1;
        }

        peer = rank - p2r;
        if(rank >= p2r && rank < pow(2, r + 1)) {
          parent = peer;
        }
      }
	//printf("rank = %d parent = %d no_childs = %d\n", rank, parent, no_childs);

	for(i=0;i<RUNS;i++){
		MPI_Barrier(MPI_COMM_WORLD);

		cdone=0; count=0;
		t1 = MPI_Wtime();
		
		// set up all recv from left and right tree
		if(rank!=root) {// if not root setup all recvs			
			MPI_Recv(msg,SIZE,MPI_CHAR,parent,0,MPI_COMM_WORLD,&stt);
			for(int l = 0; l < no_childs; l++) {
				int to_send = childs[l]; 
				MPI_Isend(msg,SIZE,MPI_CHAR,to_send,0,MPI_COMM_WORLD,&sreq[count++]);
			}	
		} 
		else {  // if root then setup all send's
			for(int l = 0; l < no_childs; l++) {
				int to_send = childs[l]; 
				MPI_Isend(outmsg,SIZE,MPI_CHAR,to_send,0,MPI_COMM_WORLD,&sreq[count++]);
			}	
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
