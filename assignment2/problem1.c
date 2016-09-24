#include <stdio.h>
#include <omp.h>

/*Calculate inner product of two vectors*/
int main () {
    int const N = 100;
    int i, k;
    double a[N], b[N]; /* declare vectors a and b */
    double dot_prod = 0.0;

    for (i=0; i<N; i++) {
	a[i] = 3.14;
	b[i] = 6.67;
    }
#pragma omp parallel
    {
#pragma omp for reduction(+: dot_prod)
	for (i = 0; i<N; i++)
	    dot_prod = dot_prod + a[i] * b[i];
    }
    printf("Inner product of a[] and b[] = %g\n", dot_prod);

    return 0;
}
