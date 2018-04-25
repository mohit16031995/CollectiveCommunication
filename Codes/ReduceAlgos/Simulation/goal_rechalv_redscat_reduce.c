#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ifroot if(rank==0)

int main(int argc,char *argv[]){
	int rank, p, root = 0, index, cdone=0;
	int i, j, k;
	long int SIZEinbytes, SIZE, CSIZE, count, thissize=0, time=0, originaltime=0;
	int l=0,r=0,mid=0,ll=0,point=0;		
	int roundsdone=0, collected=0, thisround=0, originalsize=0; 
	char* ptr;

	if(argc!=3){	 
		printf("Usage: <program> <number of proceses> <message_size> \n");
		exit(0);
	}	
	
	p = strtol(argv[1], &ptr, 10);
	SIZEinbytes = strtol(argv[2], &ptr, 10);

	//no of steps
	int rounds=0;
	rounds = (int)log2(p);
	int originalp = p;
	p = 1<<rounds;
	int extra = originalp-p;


	SIZE = SIZEinbytes/sizeof(int);		//no of elements in int array
	if (SIZE%p!=0)
		SIZE = SIZE + p - SIZE%p;
	SIZEinbytes = SIZE*sizeof(int);
	originalsize = SIZEinbytes;


	///////////////////////////////////////////////////milliseconds
	time = 3*SIZE*100;
	originaltime=time;
	///////////////////////////////////////////////////

	//2 power ith nodes on left and right
	int right[rounds];
	int left[rounds];
	int partner[rounds];
	int recvatindex[rounds];	//data to be received in a particular round is data belonging to nodes a,b,c,d for eg. then recvatindex=a
	

	printf("num_ranks %d\n", p);
	
	for (rank=0;rank<originalp;rank++)
	{
		count=1;
		time = originaltime;
		SIZE = originalsize;
		printf("\nrank %d {\n",rank);
		if (extra != 0)
		{
			if (rank >= p)
			{
				printf("l%ld: send %ldb to %d tag 0\n", count++, SIZE, rank-extra);
				printf("}\n");
				continue;
			}
			else if (rank >= p-extra)
			{
				printf("l%ld: recv %ldb from %d tag 0\n", count++, SIZE, rank+extra);
				printf("l%ld: calc %ld\n", count++, time);
				printf("l%ld requires l%ld\n", count-1, count-2);
			}
		}
		for (i=0; i<rounds; i++)
		{
			right[i] = 0;
			left[i] = 0;
		}
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
			}
			else
			{
				partner[i] = right[i];
			}
		}

		thissize = SIZE/2;

			printf("l%ld: send %ldb to %d tag 0\n", count++, thissize, partner[rounds-1]);
			if (rank >= p-extra)
			{
				printf("l%ld requires l%ld\n", count-1, count-2);
			}
			printf("l%ld: recv %ldb from %d tag 0\n", count++, thissize, partner[rounds-1]);
			printf("l%ld: calc %ld\n", count++, time);
			printf("l%ld requires l%ld\n", count-1, count-2);
			thissize/=2;
			time/=2;

		for (j=rounds-2; j>=0; j--,thissize/=2,time/=2)
		{
			//MPI_Irecv (recvmsg+thissize,thissize,MPI_INT,partner[j],j,MPI_COMM_WORLD,&req[j]);
			printf("l%ld: send %ldb to %d tag 0\n", count++, thissize, partner[j]);
			printf("l%ld requires l%ld\n", count-1, count-2);
			printf("l%ld: recv %ldb from %d tag 0\n", count++, thissize, partner[j]);
			printf("l%ld: calc %ld\n", count++, time);
			printf("l%ld requires l%ld\n", count-1, count-2);
					
		}		
		
		SIZE = thissize*2;

		int recvcount=0, sendsize=0;
		j = 0;
		while (j<rounds && partner[j]>rank)
		{
			printf("l%ld: recv %ldb from %d tag 0\n", count++, SIZE<<j, partner[j]);
			printf("l%ld requires l%ld\n", count-1, count-2);
			j++;
		}
		
		ifroot ; else 
		printf("l%ld: send %ldb to %d tag 0\n", count++, SIZE<<j, partner[j]);
			
		printf("}\n");
	}
	return 0;
}
