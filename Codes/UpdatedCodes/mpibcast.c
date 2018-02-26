#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define RUNS 100

int main(int argc,char *argv[]){
  int rank,nP, i;
  long int count;

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
  
 	long int SIZE = strtol(argv[1], &ptr, 10);
	char *outmsg = malloc(SIZE * (sizeof(char)));
	if (rank == 0) {
  		for(i=0;i<SIZE;i++)
    		outmsg[i] = 'A';
	}
  
  //for(long int s=1;s<=1024*1024;s=s*1024) {
  //if(rank==0){
  //  printf("bcast %d %ld\n",nP,s);
  //}
    for(i=0;i<RUNS;i++) {      
      MPI_Barrier(MPI_COMM_WORLD);
      t1 = MPI_Wtime();
      MPI_Bcast(outmsg, SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
      t2 = MPI_Wtime() - t1; 
	//printf("rank = %d inmsg = ", rank);
	//for (int ll = 0; ll < SIZE; ll++) {
	//	printf ("%c", inmsg[ll]);
	//}
	//printf(" outmsg = ");
	//for (int ll = 0; ll < SIZE; ll++) {
	//	printf ("%c", outmsg[ll]);
	//}
	//printf("\n");
      MPI_Reduce(&t2, &res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if(rank==0){
        printf("Run %d time %1.9lf\n", i+1,res);
      }
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    //}
  MPI_Finalize();
}
