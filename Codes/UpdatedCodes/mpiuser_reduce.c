#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define RUNS 100
void MAX_DataSet(int *in, int *inout, int *len, MPI_Datatype *datatype);
int main(int argc,char *argv[]){
  int rank,nP, i;
  long int count,s;
  
  
  char *ptr;
  MPI_Status req;
  double t1,t2,res;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nP);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Op myOp;
   MPI_Op_create((MPI_User_function*)MAX_DataSet, 1, &myOp);
	if(argc!=2){
    	if(rank==0) printf("Usage: <program> <message_size>\n");
    	exit(0);
 	 }
  
  long int SIZE = strtol(argv[1], &ptr, 10);
 	count = SIZE / (sizeof(int));
 	int *selfmsg = malloc(count * (sizeof(int)));
	int *resmsg = malloc(count * (sizeof(int)));

	
  //for(long int s=1;s<=1024*1024;s=s*1024) {
  //if(rank==0){
  //  printf("bcast %d %ld\n",nP,s);
  //}
    for(i=0;i<RUNS;i++) {      
      MPI_Barrier(MPI_COMM_WORLD);
		for (int ll=0;ll<count;ll++) {
			selfmsg[ll] = (ll + rank);
			//if(rank==0) msg[i] = selfmsg[i];	
		}
      t1 = MPI_Wtime();
      MPI_Reduce(selfmsg, resmsg, count, MPI_INT, myOp, 0, MPI_COMM_WORLD);
      t2 = MPI_Wtime() - t1; 
      MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if(rank==0){
        printf("Run %d time %1.9lf\n", i+1,res);
      }
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    //}
  MPI_Finalize();
}

void MAX_DataSet(int *in, int *inout, int *len, MPI_Datatype *datatype)
{
    for (int i = 0; i < *len; i++)
    {
        *inout = *in + *inout;
        in++;
        inout++;
    }
}
