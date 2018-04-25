#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#define RUNS 1
//#define SIZE 6
//#define CHUNK 2
//#define CSIZE SIZE/CHUNK

// Macros used in reduce collective

int main(int argc,char *argv[]){
	int rank, p, root = 0, index, cdone=0;
	int count, i, j, k, SIZEinbytes, SIZE, CSIZE, logical_chunk_no;
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
	//	rightpeers[0] = p-1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(argc!=2){	 
		if(rank==0) printf("Usage: <program> <message_size> \n");
		exit(0);
	}	

	SIZEinbytes = strtol(argv[1], &ptr, 10);
	//CHUNK = atoi(argv[2]);
	//CSIZE = SIZE/(CHUNK*(sizeof(int))); 
	SIZE = SIZEinbytes/sizeof(int);		//no of elements in int array
	//printf("SIZE: %d\n",SIZE);
	int *start = malloc(SIZE*2 * sizeof(int));
	int *data[2];
	//for (i=0; i<; i++)
	//{
	//	data[i] = start+(i*SIZE);
	//}
	data[0] = start;
	data[1] = start + SIZE;
	int *selfmsg;
	selfmsg = data[0];
//	int ready[CHUNK];
	
	//no of steps
	int rounds=0;
	for (i=0,j=p; j>1; j/=2,i++);	//no of rounds. log base 2 p.
	rounds = i;
	//printf("i = %d",i);

	// for  recvs 
	MPI_Status stt;		//no of send/recv = no of rounds
	MPI_Request req[rounds];
	for (i=0;i<rounds;i++)
		req[i] = MPI_REQUEST_NULL;
	// for  send
	MPI_Status sstt[rounds];
	MPI_Request sreq[rounds];
	for (i=0;i<rounds;i++)
		sreq[i] = MPI_REQUEST_NULL;

	double t1,t2,res;

	//2 power ith nodes on left and right
	int right[rounds];
	int left[rounds];
	int partner[rounds];
	int recvatindex[rounds];	//data to be received in a particular round is data belonging to nodes a,b,c,d for eg. then recvatindex=a
	right[0]=(rank+1<p) ? rank+1 : rank+1-p;
	left[0]=(rank-1>=0) ? rank-1 : rank-1+p;
	int offset=1, prev=right[0];
	for (i=1; i<rounds; i++)
	{
		right[i] = prev+offset;
		offset*=2;
		if (right[i]>p-1)
			right[i]-=p;
		prev = right[i];
	}
	offset=1, prev=left[0];
	for (i=1; i<rounds; i++)
	{
		left[i] = prev-offset;
		offset*=2;
		if (left[i]<0)
			left[i]+=p;
		prev = left[i];
	}
	int r=rank;
	for (i=0; i<rounds; i++,r=r>>1)
	{
		if (r&1) //odd
		{
			partner[i] = left[i];
			recvatindex[i] = (r-1)<<i;
		}
		else
		{
			partner[i] = right[i];
			recvatindex[i] = (r+1)<<i;
		}
	}

	for (i=0;i<RUNS;i++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		count=0;
		int roundsdone=0, collected = rank, thisround=0; //collected represents location of the data we have collected till now
		for (int ll=0; ll<SIZE*p; ll++)
			start[i]=0;
		//printf("Selfmsg of rank %d\n",rank);
		for (int ll=0;ll<SIZE;ll++) {
			selfmsg[ll] = rank;
			//printf("%d ",*(data[0]+ll));
			//if(rank==0) msg[i] = selfmsg[i];	
		}
		//printf("\n\n");
		for (j=0; j<rounds; j++)
		{
			MPI_Irecv (data[1],SIZE,MPI_INT,partner[j],j,MPI_COMM_WORLD,&req[j]);
			printf("rank %d to recv %d elements from %d at index %d\n",rank,SIZE,partner[j],SIZE);
		}
		MPI_Isend (data[0],SIZE,MPI_INT,partner[0],0,MPI_COMM_WORLD,&sreq[count++]);
		roundsdone=1;
		while (roundsdone<rounds)
		{
			MPI_Wait(&req[roundsdone-1], &stt);
			//printf("ok waiting, rank %d",rank);
		/*	if(thisround == MPI_UNDEFINED) 
			{
				printf("Unexpected error!\n");
				MPI_Abort(MPI_COMM_WORLD, 1);
			}*/
			thisround=roundsdone-1;
			//if (collected>recvatindex[thisround]) 
			//	collected = recvatindex[thisround];
			printf ("rank %d received round %d\t collected= %d\n",rank,thisround,collected);
			for(int i = 0; i < SIZE; i++)
				*(data[0]+i) = *(data[0]+i) + *(data[1]+i);
			if (rank==0)
			{
				//printf("Rank : %d\n\n",rank);
				int ii;
				for (ii=0; ii<SIZE*2; ii++)
					printf ("%d  ",start[ii]);
				printf("\n\n");
			}
			MPI_Isend (data[0],SIZE,MPI_INT,partner[thisround+1],thisround+1,MPI_COMM_WORLD,&sreq[count++]);
			printf("rank %d to send %d elements to %d at index %d\n",rank,SIZE,partner[thisround+1],collected);
			roundsdone++;
		}

		//for(int i = 0; i < SIZE; i++)
		//	*(data[0]+i) = *(data[0]+i) + *(data[1]+i);
		
		//MPI_Waitall(rounds,req,sstt);
		MPI_Wait(&req[roundsdone-1], &stt);

		for(int i = 0; i < SIZE; i++)
			*(data[0]+i) = *(data[0]+i) + *(data[1]+i);

		MPI_Waitall(count,sreq,sstt);  // wait for  all send to finish
		
////		printf("My work is done here, rank %d\n",rank);
		t2 = MPI_Wtime() - t1;
		MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

		//		for (i=0; i<SIZE; i++)
		//			printf ("%d  ",selfmsg[i]);

		if(rank==0)
		{
			printf("Run %ld time %1.9lf\n", i+1,res);
			for (int ii=0; ii<SIZE*2; ii++)
				printf ("%d  ",start[ii]);
			printf("\n");
		}

	}
	MPI_Finalize();
}
