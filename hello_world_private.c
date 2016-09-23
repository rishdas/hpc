#include <stdio.h>
#include <omp.h>

int main() {
    int nthreads, tid;
#pragma omp parallel private(nthreads, tid)
    {
	tid = omp_get_thread_num();
	printf("Hello World from thread no: %d\n", tid);
	if (tid == 0) {
	    nthreads = omp_get_num_threads();
	    printf("Total no of threads: %d\n", nthreads);
	}
    }
}
