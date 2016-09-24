#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[])
{
    const int N = 1000;
    int chunk = 10;
    double mat[N][N], vect[N], result[N];
    int i, j;
    int sum = 0;
    //Init Matrix
    for (i=0; i<N; i++)
	for (j = 0; j<N; j++)
	    mat[i][j] = 1.3*i*j;

    //Init Vector
    for (i = 0; i<N; i++)
	vect[i] = 2.7*i;

    //Init result
    for (i = 0; i<N; i++)
	result[i] = 0.0;

#pragma omp parallel shared(sum)
    {
#pragma omp for schedule(dynamic, chunk)
	for (i=0; i<N; i++)
	    for (j = 0; j<N; j++)
		result[i] = result[i] + mat[i][j]*vect[j];
    }
    for (i = 0; i<N; i++)
	printf("%g\n", result[i]);
}
