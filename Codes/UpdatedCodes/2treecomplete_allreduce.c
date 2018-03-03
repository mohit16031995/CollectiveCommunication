#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#define RUNS 100
//#define SIZE 6
//#define CHUNK 2
//#define CSIZE SIZE/CHUNK

// Macros used in reduce collective

int main(int argc,char *argv[]){
	int rank, p, root = 0, index, cdone=0;
	long int count, i, j, k, SIZE, CSIZE, logical_chunk_no;
	long int CHUNK;
	char* ptr;

	int parentLeft = -1;
	int parentRight = -1;

	int leftChildren = 0;
	int rightChildren = 0;

	int leftPeers[2];
	int rightPeers[2];


	int seed = time(NULL);
	srand(seed);
	//char inmsg[SIZE], selfmsg[SIZE];

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	//	righteers[0] = p-1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(argc!=3){	 
		if(rank==0) printf("Usage: <program> <message_size> <nchunk>\n");
		exit(0);
	}	

	SIZE = strtol(argv[1], &ptr, 10);
	CHUNK = atoi(argv[2]);
	CSIZE = SIZE/(CHUNK*(sizeof(int))); 
	SIZE = CSIZE*CHUNK;

	int *msg1 = malloc((SIZE*2+1) * sizeof(int));
	int *selfmsg = malloc((SIZE+1) * sizeof(int));
	int *msg2 = msg1+SIZE;
	int ready[CHUNK];
	int *Reducedmsg = malloc((SIZE+1) * sizeof(int));
	// for  recvs 
	MPI_Status stt;
	MPI_Request *req1= calloc(CHUNK*2,sizeof(MPI_Request));
	MPI_Request *req2 = req1+CHUNK;
	//printf("I am rank %d, my req1 is at %p\n",rank,req1);
	for (i=0;i<CHUNK*2;i++)
		req1[i] = MPI_REQUEST_NULL;
	// for  send
	MPI_Status sstt[CHUNK*2];
	MPI_Request sreq[CHUNK*2];
	MPI_Request req[CHUNK];

	for (i=0;i<CHUNK;i++)
		sreq[i] = MPI_REQUEST_NULL;

	double t1,t2,res;

	for (i=0;i<SIZE;i++) {
		selfmsg[i] = (i%100 + rank%100)%100;
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
	else
	{
		leftChildren = 1;
		rightChildren = 1;
		leftPeers[0] = 1;
		rightPeers[0] = p-1;
	}

	//	double timings[2][50][515];

	for (i=0;i<RUNS;i++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		for (int ll=0;ll<SIZE;ll++) {
			selfmsg[ll] = (ll%100 + rank%100)%100;
			//if(rank==0) msg[i] = selfmsg[i];	
		}
		cdone=0; count=0;
		for (k=0;k<CHUNK;k++)
			ready[k]=0;
		t1 = MPI_Wtime();

		// set up all recv from left and right tree
		if (leftChildren)				//if not leafleft setup all even recvs
			for (j=0;j<CHUNK;j+=2) 
			{
				MPI_Irecv(msg1+j*CSIZE,CSIZE,MPI_INT,leftPeers[0],j,MPI_COMM_WORLD,&req1[j]);
				if (leftChildren==2)
					MPI_Irecv(msg2+j*CSIZE,CSIZE,MPI_INT,leftPeers[1],j,MPI_COMM_WORLD,&req1[j+CHUNK]);
				//if (j+CHUNK==CHUNK || j+CHUNK==CHUNK+1) 
				//printf("\tI am rank %d, recv tag j %d from left tree child %d\n", rank,j,leftPeers[1]);
			}
		if (rightChildren)				//if not leafRight setup all odd recvs
			for (j=1;j<CHUNK;j+=2) 
			{
				MPI_Irecv(msg1+j*CSIZE,CSIZE,MPI_INT,rightPeers[0],j,MPI_COMM_WORLD,&req1[j]);
				if (rightChildren==2)
					MPI_Irecv(msg2+j*CSIZE,CSIZE,MPI_INT,rightPeers[1],j,MPI_COMM_WORLD,&req1[j+CHUNK]);
				//if (j+CHUNK==CHUNK || j+CHUNK==CHUNK+1) 
				//printf("\tI am rank %d, recv tag j %d from right tree child %d\n", rank,j,rightPeers[1]);
			}

		// setup all sends in left and right leaves
		if (!leftChildren || !rightChildren)
			for (j=0;j<CHUNK;j++) 
			{	
				if (!leftChildren && !(j%2)) 
					MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,parentLeft,j,MPI_COMM_WORLD,&sreq[count++]);
				if (!rightChildren && (j%2))
					MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,parentRight,j,MPI_COMM_WORLD,&sreq[count++]);
			}

		//MPI_Waitall(count,sreq,sstt);  // wait for  all send to finish
		int sent_as_leaf = count;                 
////		printf ("%d leaf role completed\n",rank);
		// check which recv finish and setup send for them
		int chunks_to_recv = ((CHUNK+1)/2)*leftChildren + (CHUNK/2)*rightChildren;
		//chunks_to_recv = CHUNK;
		//printf ("%d, left=%d right=%d, rank %d\n",chunks_to_recv,leftChildren,rightChildren,rank);
		if(rank!=root) 
		{
			while(cdone < (chunks_to_recv)) 
			{		
				MPI_Waitany(CHUNK*2, req1, &index, &stt);
				//printf("ok waiting, rank %d",rank);
				if(index == MPI_UNDEFINED) 
				{
					printf("Unexpected error!\n");
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				logical_chunk_no = (index>=CHUNK) ? (index-CHUNK) : index;
////				printf ("\trank %d received chunk %d (actual %d)\n",rank,logical_chunk_no,index);
				
				double t11 = MPI_Wtime();
				j=index*CSIZE;
				
					for (k=logical_chunk_no*CSIZE;k<(logical_chunk_no+1)*CSIZE;k++) 
					{
						selfmsg[k] = (msg1[j]+selfmsg[k]);
						j++;
					}
				double t12 = MPI_Wtime()-t11;
////				printf("Calculation time 1 %d %1.9f\n", CSIZE, t12);
				ready[logical_chunk_no]++;
				if ( (logical_chunk_no%2) ? ready[logical_chunk_no]==rightChildren : ready[logical_chunk_no]==leftChildren )
				{	//forward this chunk ;
////					printf("chunk no. %d is ready to forward from rank %d\n",logical_chunk_no,rank);
					
					if (!(logical_chunk_no%2))
						MPI_Isend(selfmsg+logical_chunk_no,CSIZE,MPI_INT,parentLeft,logical_chunk_no,MPI_COMM_WORLD,&sreq[count++]);
					else
						MPI_Isend(selfmsg+logical_chunk_no,CSIZE,MPI_INT,parentRight,logical_chunk_no,MPI_COMM_WORLD,&sreq[count++]);
				}
				cdone++;
				//printf("cdone = %d, rank %d",cdone,rank);
			}
////			printf("done forwarding %d\n",rank);
		}
		else
			while(cdone < (chunks_to_recv)) 
			{		
				MPI_Waitany(CHUNK*2, req1, &index, &stt);
				if(index == MPI_UNDEFINED) 
				{
					printf("Unexpected error!\n");
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				logical_chunk_no = (index>=CHUNK) ? (index-CHUNK)*CSIZE : index*CSIZE;
				double t11 = MPI_Wtime();
				j=index*CSIZE;
				
					for (k=logical_chunk_no;k<logical_chunk_no+CSIZE;k++) 
					{
						selfmsg[k] = (msg1[j]+selfmsg[k]);		
						j++;
					}
			
				double t12 = MPI_Wtime() - t11;
////				printf("Calculation time 2 %d %1.9f\n", CSIZE, t12);	
				cdone++;
			}



		MPI_Waitall(count,sreq,sstt);  // wait for  all send to finish
////		printf("My work is done here, rank %d\n",rank);
		if (rank == 0) {
			//printf("Reduce Completed\n");
		}
				///////////////////////////////////////////////////////////////////////////BROADCAST////////////////////////////////////////
		cdone=0; count=0;
		// set up all recv from left and right tree
		if(rank!=root) {// if not root setup all recvs
		   for(j=0;j<CHUNK;j++) {
			// left tree				
			if(!(j%2)) MPI_Irecv(Reducedmsg+j*CSIZE,CSIZE,MPI_INT,parentLeft,CHUNK+j,MPI_COMM_WORLD,&req[j]);			
			// right tree
			else MPI_Irecv(Reducedmsg+j*CSIZE,CSIZE,MPI_INT,parentRight,CHUNK+j,MPI_COMM_WORLD,&req[j]);
		   }
		} else {  // if root then setup all send's
		   for(j=0;j<CHUNK;j++) {
			// left tree				
			if(!(j%2)) { 
				MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,1,CHUNK+j,MPI_COMM_WORLD,&sreq[count++]);
			}			
			// right tree
			else { 
				MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,p-1,CHUNK+j,MPI_COMM_WORLD,&sreq[count++]);
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
		//printf("22rank = %d chunk recvd = %d\n", rank, index);
			index = index;
			if(index%2) {  // right tree
				if(rightChildren) {
					MPI_Isend(Reducedmsg+index*CSIZE,CSIZE,MPI_INT,rightPeers[0],CHUNK+index,MPI_COMM_WORLD,&sreq[count++]);
				}
				if(rightChildren==2) {
					MPI_Isend(Reducedmsg+index*CSIZE,CSIZE,MPI_INT,rightPeers[1],CHUNK+index,MPI_COMM_WORLD,&sreq[count++]);
				}
			} else {
				if(leftChildren) {
					MPI_Isend(Reducedmsg+index*CSIZE,CSIZE,MPI_INT,leftPeers[0],CHUNK+index,MPI_COMM_WORLD,&sreq[count++]);
				}
				if(leftChildren==2) {
					MPI_Isend(Reducedmsg+index*CSIZE,CSIZE,MPI_INT,leftPeers[1],CHUNK+index,MPI_COMM_WORLD,&sreq[count++]);
				}
			}
			cdone++;
		  }

		//printf("rank = %d\n", rank);
		MPI_Waitall(count,sreq,sstt);  // wait for all send to finish
			
		t2 = MPI_Wtime() - t1;
		MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

		//		for (i=0; i<SIZE; i++)
		//			printf ("%d  ",selfmsg[i]);

		if(rank==0)
		{
			printf("Run %ld time %1.9lf\n", i+1,res);
		//	for (int ii=0; ii<SIZE; ii++) {
		//		printf("%d ", selfmsg[ii]);
		//	}
		//	printf("\n");
		}

	}
	MPI_Finalize();
} 

