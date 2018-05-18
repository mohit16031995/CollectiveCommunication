#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#define RUNS 1
#define ifroot if(rank==0)
//#define SIZE 6
//#define CHUNK 2
//#define CSIZE SIZE/CHUNK

// Macros used in reduce collective

int main(int argc,char *argv[]){
	int rank, p, root = 0, index, cdone=0;
	int count, i, j, k, SIZEinbytes, SIZE, CSIZE, logical_chunk_no;
	int l=0,r=0,mid=0,recvsize=0,ll=0,point=0;		
	int roundsdone=0, collected=0, thisround=0; 
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

	int selfmsg[SIZE];
	int recvmsg[SIZE];
//	int ready[CHUNK];
	
	//no of steps
	int rounds=0;
	for (i=0,j=p; j>1; j=j>>1,i++);	//no of rounds. log base 2 p.
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

	double t1,t2,res,time=0;

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
	int rr=rank;
	for (i=0; i<rounds; i++,rr=rr>>1)
	{
		if (rr&1) //odd
		{
			partner[i] = left[i];
			recvatindex[i] = (rr-1)<<i;
		}
		else
		{
			partner[i] = right[i];
			recvatindex[i] = (rr+1)<<i;
		}
	}

	for (i=0;i<RUNS;i++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		count=0;
		l=0,r=SIZE-1,mid=(l+r+1)/2,point=0;		//left to mid-1 ,, mid to right ,, mid=(left+right+1)/2
		roundsdone=0, collected = rank, thisround=0; //collected represents location of the data we have collected till now
		for (int ii=0;ii<SIZE;ii++) {
			selfmsg[ii] = rank;
		}
		t1 = MPI_Wtime();
		recvsize = SIZE>>1;
		for (j=rounds-1; j>=0; j--,recvsize=recvsize>>1)
		{
			MPI_Irecv (recvmsg+recvsize,recvsize,MPI_INT,partner[j],j,MPI_COMM_WORLD,&req[j]);
			//printf("rank %d to recv %d elements from %d at index %d at tag %d\n",rank,SIZE<<(rounds-j-1),partner[j],recvatindex[j],j);
		}
		if (partner[rounds-1]<rank) 
		{	
			ll=mid;
			point=l;
			l=mid;
			mid=(l+r+1)/2;
		}
		if (partner[rounds-1]>rank) 
		{
			ll=l;
			point=mid;
			r=mid-1;
			mid=(l+r+1)/2;
		}
ifroot printf("ll=%d",ll);
			
		MPI_Isend (selfmsg+point,SIZE>>1,MPI_INT,partner[rounds-1],rounds-1,MPI_COMM_WORLD,&sreq[count++]);
//		printf("rank %d to send %d elements to %d at index %d at tag %d\n",rank,SIZE,partner[rounds-1],collected,rounds-1);		
		roundsdone=1;
		recvsize = SIZE>>1;
		for (int torecv=rounds-1; torecv>0; torecv--,recvsize=recvsize>>1)
		{
			MPI_Wait(&req[torecv], &stt);
			for (int k=0,t=ll; k<recvsize; k++,t++)
			{
				selfmsg[t]+=recvmsg[recvsize+k];
			}
			ifroot
				printf("\nll=%d rr=%d\n",ll,ll+recvsize);
			if (partner[torecv-1]<rank) 
			{	
				ll=mid;
				point=l;
				l=mid;
				mid=(l+r+1)/2;
			}
			else
			{
				ll=l;
				point=mid;
				r=mid-1;
				mid=(l+r+1)/2;
			}
			ifroot
			{
				for (int ii=0; ii<SIZE; ii++)
					printf("%d ",selfmsg[ii]);
				printf("Recvd round %d\n",torecv);
			}
			MPI_Isend (selfmsg+point,recvsize>>1,MPI_INT,partner[torecv-1],torecv-1,MPI_COMM_WORLD,&sreq[count++]);
		}		
		MPI_Waitall(rounds,req,sstt);
	
		for (int k=0,t=ll; k<recvsize; k++,t++)
		{
			selfmsg[t]+=recvmsg[recvsize+k];
		}
	
		MPI_Waitall(count,sreq,sstt);  // wait for  all send to finish
		
		t2 = MPI_Wtime() - t1;
		MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		time+=res;
				
		printf("rank %d run %d total %d\n",rank,i+1,selfmsg[ll]);

		if(rank==0)
		{
			printf("Run %ld time %1.9lf\n", i+1,res);
		//	for (int ii=0; ii<SIZE*p; ii++)
		//		printf ("%d  ",start[ii]);
		//	printf("\n");
		}

	}
	if (rank==0) 
		printf("Avg time = %lf",time/RUNS);
	MPI_Finalize();
} 

