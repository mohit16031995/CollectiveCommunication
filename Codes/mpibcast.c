#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define RUNS 100
#define SIZE 1048576
int main(int argc,char *argv[]){
  int rank,nP;
  long int count,i,s;
  char inmsg[SIZE], outmsg[SIZE];
  char *ptr;
  MPI_Status req;
  double t1,t2,res;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nP);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(argc!=2){
    if(rank==0) printf("Usage: <program> <message_size>\n");
    exit(0);
  }
  
  s = strtol(argv[1], &ptr, 10);

  for(i=0;i<SIZE;i++)
    outmsg[i] = 'A';
  
  //for(long int s=1;s<=1024*1024;s=s*1024) {
  //if(rank==0){
  //  printf("bcast %d %ld\n",nP,s);
  //}
    for(i=0;i<RUNS;i++) {      
      MPI_Barrier(MPI_COMM_WORLD);
      t1 = MPI_Wtime();
      MPI_Bcast(outmsg, s, MPI_CHAR, 0, MPI_COMM_WORLD);
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
