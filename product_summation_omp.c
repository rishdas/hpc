#include <stdio.h>
#include <omp.h>

int main() {
    const int N = 1000000;
    int i;
    double a[N], b[N], dot_prod;

    dot_prod = 0.0;
    for (i = 0; i<N; i++) {
	a[i] = 6.67;
	b[i] = 3.14;
    }
#pragma omp parallel
    {
#pragma omp for reduction(+: dot_prod)
	for (i = 0; i<N; i++)
	    dot_prod = dot_prod + (a[i] * b[i]);
    }
    printf(" Dot product = %g\n", dot_prod);

    return 0;
}
