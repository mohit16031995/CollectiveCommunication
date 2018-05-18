#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#define SIZE 1024*1000

void main() {
	int *a=(void*)malloc(SIZE*sizeof(int));
	int *b=(void*)malloc(SIZE*sizeof(int));
//	float *c=(void*)malloc(SIZE*sizeof(float));
	struct timeval  tv1, tv2;

	/* single thread */
	gettimeofday(&tv1, NULL);	
 	for(int n=0; n<SIZE; ++n)
 	{
	   a[n] = a[n]+b[n];
	}
	gettimeofday(&tv2, NULL);

	printf ("Total time one thread = %lf seconds\n",
         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double) (tv2.tv_sec - tv1.tv_sec));

	gettimeofday(&tv1, NULL);
	for (int n = 0; n < SIZE; ++n) {
		a[n] = (a[n]%100+b[n]%100)%100;
	}	
	gettimeofday(&tv2, NULL);	
	printf ("Total time with modulo = %lf seconds\n",
         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double) (tv2.tv_sec - tv1.tv_sec));
	/* multiple threads */
	gettimeofday(&tv1, NULL);
	#pragma omp for
 	for(int n=0; n<SIZE; ++n)
 	{
	   a[n] = a[n]+b[n];
	}
	gettimeofday(&tv2, NULL);

	printf ("Total time = %lf seconds\n",
         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
         (double) (tv2.tv_sec - tv1.tv_sec));
}


