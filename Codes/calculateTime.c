#include "mpi.h"
#include "mpi.h"
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

int main(int argc,char *argv[]){
	int seed = time(NULL);
	srand(seed);
	//char inmsg[SIZE], selfmsg[SIZE];
	int p, rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	//	righteers[0] = p-1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	for (int CSIZE = 1024; CSIZE < 104857; CSIZE+=1024){
		int *msg1 = malloc(CSIZE*2 * sizeof(int));
		int *selfmsg = malloc(CSIZE * sizeof(int));

		for (int i=0;i<CSIZE;i++) {
			selfmsg[i] = (i%100 + rank%100)%100;
			//if(rank==0) msg[i] = selfmsg[i];	
		}
		double avg = 0;
		int Runs = 10;
		for (int l = 0; l < Runs; l++) {
			double t1 = MPI_Wtime();
			#pragma omp parallel
			{
						#pragma omp for
						for (int k = 0; k < CSIZE; k++) {
							selfmsg[k] = (msg1[k]%100+selfmsg[k]%100)%100;
						}
			}
			double t2 = MPI_Wtime()-t1;
			avg += t2;
			printf("Run %d rank %d CSIZE %d time %1.9f\n", l+1, rank, CSIZE, t2);
		}
		avg = avg / Runs;
		printf ("rank %d CSIZE %d Avgtime %1.9f\n", rank, CSIZE, avg);
	}
	MPI_Finalize();
	return 0;
}
