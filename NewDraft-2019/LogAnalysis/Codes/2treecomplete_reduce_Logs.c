#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#define RUNS 20
//#define SIZE 6
//#define CHUNK 2
//#define CSIZE SIZE/CHUNK
#define toLog(i) (i!=0 && i%(RUNS/4) == 0)
//#define toLog(i) (i==0)
#define logIRecv(t,chunk_no,from) printf("LogIRecv: %1.9f rank %d called irecv for chunk %d from %d in run %ld %d_%ld_%d\n",t,rank,chunk_no,from,i+1,p,sizeInBytes,CHUNK);
#define logISend(t,chunk_no,to) printf("LogISend: %1.9f rank %d called isend for chunk %d to %d in run %ld %d_%ld_%d\n",t,rank,chunk_no,to,i+1,p,sizeInBytes,CHUNK);
#define logEndTime(t) printf("LogEndTime: %1.9f rank %d finished run %ld %d_%ld_%d\n",t,rank,i+1,p,sizeInBytes,CHUNK);
#define logReceived(t,logical_chunk_no,index,from) printf("LogReceived: %1.9f rank %d received chunk %d (%d) from %d in run %ld %d_%ld_%d\n",t,rank,logical_chunk_no,index,from,i+1,p,sizeInBytes,CHUNK);

// Macros used in reduce collective

int main(int argc,char *argv[]){
	int rank, p, root = 0, index, cdone=0;
	long int count, i, j, k, SIZE, CSIZE, logical_chunk_no;
	int CHUNK;
	char* ptr;

	int parentLeft = -1;
	int parentRight = -1;

	int leftChildren = 0;
	int rightChildren = 0;

	int leftPeers[2];
	int rightPeers[2];

	leftPeers[0] = -1;
	leftPeers[1] = -1;

	rightPeers[0] = -1;
	rightPeers[1] = -1;


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
	CSIZE = ceil(SIZE/(CHUNK*(sizeof(int))*1.0)); 
	long int sizeInBytes = SIZE;
	SIZE = CSIZE*CHUNK;

	int *msg1 = malloc(SIZE*2 * sizeof(int));
	int *selfmsg = malloc(SIZE * sizeof(int));
	int *msg2 = msg1+SIZE;
	int ready[CHUNK];

	// for  recvs 
	MPI_Status stt;
	MPI_Request *req1= calloc(CHUNK*2,sizeof(MPI_Request));
	MPI_Request *req2 = req1+CHUNK;
	//printf("I am rank %d, my req1 is at %p\n",rank,req1);
	for (i=0;i<CHUNK*2;i++)
		req1[i] = MPI_REQUEST_NULL;
	// for  send
	MPI_Status sstt[CHUNK];
	MPI_Request sreq[CHUNK];
	for (i=0;i<CHUNK;i++)
		sreq[i] = MPI_REQUEST_NULL;

	double t1,t2,res;

	
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

	
	printf("Tree: rank = %d, leftParent = %d, rightParent = %d, leftleftChild = %d, leftRightChild = %d, rightLeftChild = %d, rightRightChild = %d\n", rank, parentLeft, parentRight, leftPeers[0], leftPeers[1], rightPeers[0], rightPeers[1]);
	
	
	long int chunks_to_recv = ((CHUNK+1)/2)*leftChildren + (CHUNK/2)*rightChildren;

	//for logging
	double recvd_timestamps[chunks_to_recv];
	double irecv_timestamps[CHUNK][2];	//0 for left child, 1 for right child
	double isend_timestamps[CHUNK];
	int logicalNos[chunks_to_recv];
	int indexes[chunks_to_recv];
	
	
	for (i=0;i<RUNS;i++)
	{
		for (int ll=0;ll<SIZE;ll++) {
			selfmsg[ll] = 1;
			//if(rank==0) msg[i] = selfmsg[i];	
		}

		//logs resetting
		for (int ii=0; ii<chunks_to_recv; ii++) {
			recvd_timestamps[ii]=0;
			logicalNos[ii]=0;
			indexes[ii]=0;
		}
		for (int ii=0; ii<CHUNK; ii++) {
			irecv_timestamps[ii][0]=0;
			irecv_timestamps[ii][1]=0;
			isend_timestamps[ii]=0;
		}
		
		cdone=0; count=0;
		for (k=0;k<CHUNK;k++)
			ready[k]=0;

		MPI_Barrier(MPI_COMM_WORLD);
		t1 = MPI_Wtime();

		// set up all recv from left and right tree
		if (leftChildren)				//if not leafleft setup all even recvs
			for (j=0;j<CHUNK;j+=2) 
			{
				MPI_Irecv(msg1+j*CSIZE,CSIZE,MPI_INT,leftPeers[0],j,MPI_COMM_WORLD,&req1[j]);
				irecv_timestamps[j][0] = MPI_Wtime();
				if (leftChildren==2){
					MPI_Irecv(msg2+j*CSIZE,CSIZE,MPI_INT,leftPeers[1],j,MPI_COMM_WORLD,&req1[j+CHUNK]);					
					irecv_timestamps[j][1] = MPI_Wtime();
				}
				//if (j+CHUNK==CHUNK || j+CHUNK==CHUNK+1) 
				//printf("\tI am rank %d, recv tag j %d from left tree child %d\n", rank,j,leftPeers[1]);
			}
		if (rightChildren)				//if not leafRight setup all odd recvs
			for (j=1;j<CHUNK;j+=2) 
			{
				MPI_Irecv(msg1+j*CSIZE,CSIZE,MPI_INT,rightPeers[0],j,MPI_COMM_WORLD,&req1[j]);				
				irecv_timestamps[j][0] = MPI_Wtime();
				if (rightChildren==2){
					MPI_Irecv(msg2+j*CSIZE,CSIZE,MPI_INT,rightPeers[1],j,MPI_COMM_WORLD,&req1[j+CHUNK]);
					irecv_timestamps[j][1] = MPI_Wtime();
				}
				//if (j+CHUNK==CHUNK || j+CHUNK==CHUNK+1) 
				//printf("\tI am rank %d, recv tag j %d from right tree child %d\n", rank,j,rightPeers[1]);
			}

		// setup all sends in left and right leaves
		if (!leftChildren || !rightChildren)
			for (j=0;j<CHUNK;j++) 
			{	
				if (!leftChildren && !(j%2)){
					MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,parentLeft,j,MPI_COMM_WORLD,&sreq[count++]);
					isend_timestamps[j] = MPI_Wtime();
				}
				if (!rightChildren && (j%2)){
					MPI_Isend(selfmsg+j*CSIZE,CSIZE,MPI_INT,parentRight,j,MPI_COMM_WORLD,&sreq[count++]);
					isend_timestamps[j] = MPI_Wtime();
				}
			}

		//MPI_Waitall(count,sreq,sstt);  // wait for  all send to finish
////		printf ("%d leaf role completed\n",rank);
		// check which recv finish and setup send for them
		//chunks_to_recv = CHUNK;
		//printf ("%d, left=%d right=%d, rank %d\n",chunks_to_recv,leftChildren,rightChildren,rank);
		cdone = 0;
		if(rank!=0) 
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
				if (toLog(i)) {
					recvd_timestamps[cdone] = MPI_Wtime();
					logicalNos[cdone] = logical_chunk_no;
					indexes[cdone] = index;
					//printf ("Log: Time %1.9lf Rank %d received chunk %ld (actual %d)\n", MPI_Wtime()-t1,rank,logical_chunk_no,index);
				}
				
				//double t11 = MPI_Wtime();
				j=index*CSIZE;
				
					for (k=logical_chunk_no*CSIZE;k<(logical_chunk_no+1)*CSIZE;k++) 
					{
						selfmsg[k] = (msg1[j]+selfmsg[k]);
						j++;
					}
				//double t12 = MPI_Wtime()-t11;
////				printf("Calculation time 1 %d %1.9f\n", CSIZE, t12);
				ready[logical_chunk_no]++;
				if ( (logical_chunk_no%2) ? ready[logical_chunk_no]==rightChildren : ready[logical_chunk_no]==leftChildren )
				{	//forward this chunk ;
////					printf("chunk no. %d is ready to forward from rank %d\n",logical_chunk_no,rank);
					
					if (!(logical_chunk_no%2)){
						MPI_Isend(selfmsg+logical_chunk_no*CSIZE,CSIZE,MPI_INT,parentLeft,logical_chunk_no,MPI_COMM_WORLD,&sreq[count++]);
						isend_timestamps[logical_chunk_no] = MPI_Wtime();
					}
					else {
						MPI_Isend(selfmsg+logical_chunk_no*CSIZE,CSIZE,MPI_INT,parentRight,logical_chunk_no,MPI_COMM_WORLD,&sreq[count++]);
						isend_timestamps[logical_chunk_no] = MPI_Wtime();
					}
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
				if (toLog(i)) {
					recvd_timestamps[cdone] = MPI_Wtime();
					logicalNos[cdone] = logical_chunk_no/CSIZE;
					indexes[cdone] = index;
					//printf ("Log: Time %1.9lf Rank %d received chunk %ld (actual %d)\n", MPI_Wtime()-t1,rank,logical_chunk_no,index);
				}
				//double t11 = MPI_Wtime();
				j=index*CSIZE;
				
					for (k=logical_chunk_no;k<logical_chunk_no+CSIZE;k++) 
					{
						selfmsg[k] = (msg1[j]+selfmsg[k]);		
						j++;
					}
			
				//double t12 = MPI_Wtime() - t11;
//				printf("Calculation time 2 %d %1.9f\n", CSIZE, t12);	
				cdone++;
			}



		MPI_Waitall(count,sreq,sstt);  // wait for  all send to finish
		
////		printf("My work is done here, rank %d\n",rank);
		t2 = MPI_Wtime() - t1;
		MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		

		if (toLog(i)){
			for (int ii=0; ii<CHUNK; ii++) {
				if (irecv_timestamps[ii][0]!=0)
					logIRecv(irecv_timestamps[ii][0]-t1,ii,ii%2==0?leftPeers[0]:rightPeers[0]);
				if (irecv_timestamps[ii][1]!=0)
					logIRecv(irecv_timestamps[ii][1]-t1,ii,ii%2==0?leftPeers[1]:rightPeers[1]);
			}
			if (rank!=0)
			for (int ii=0; ii<CHUNK; ii++) {
				logISend(isend_timestamps[ii]-t1,ii,ii%2==0?parentLeft:parentRight);
			}
			for (int ii=0,from=-1; ii<chunks_to_recv; ii++) {
				from=-1;
				if (logicalNos[ii]%2==0)
					if (logicalNos[ii]==indexes[ii])
						from = leftPeers[0];
					else
						from = leftPeers[1];
				else
					if (logicalNos[ii]==indexes[ii])
						from = rightPeers[0];
					else
						from = rightPeers[1];
				logReceived(recvd_timestamps[ii]-t1,logicalNos[ii],indexes[ii],from);
			}
			logEndTime(t2);
		}

		//Correctness Check
		if(rank==0)
		{
			printf("Run %ld time %1.9lf %d_%ld_%d\n", i+1,res,p,sizeInBytes,CHUNK);
			long int ii=0;
			for (ii=0; ii<SIZE; ii++){
				if (selfmsg[ii]!=p){
					ii = SIZE+1;
					printf("Fatal Error! Incorrect result at index %ld. Abort mission. %d_%ld_%d\n",ii,p,sizeInBytes,CHUNK);
				}
			}
			if (ii != SIZE+1) {
				printf("Result Correct run %ld %d_%ld_%d\n", i+1,p,sizeInBytes,CHUNK);
			}
		}

	}
	MPI_Finalize();
} 

