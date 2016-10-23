#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

double square(double x)
{
    return x*x;
}
double rand_M_N(double M, double N)
{
    return M + (rand() / ( RAND_MAX / (N-M) ) ) ;  
}
double try(int local_p_size)
{
    double x, y, pi;
    int hit = 0;
    for (int i = 0; i < local_p_size; i++)  {
	x = rand_M_N(-1.0, 1.0);
	y = rand_M_N(-1.0, 1.0);

	if ((square(x) + square(y)) <= 1.0)
	    hit++;
    }

    pi = 4.0 * (double)hit/(double)local_p_size;
    return pi;
}

int main (int argc, char *argv[])
{
    double	local_pi, tot_pi, pi;
    int         trys = 64000000;
    int	        rank, num_procs, ret_val, i;
    MPI_Status  status;
    int         root = 0;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srandom(rank);
    local_pi = try(trys/num_procs);

    MPI_Reduce(&local_pi, &tot_pi, 1, MPI_DOUBLE, MPI_SUM,
		    root, MPI_COMM_WORLD);
    if (rank == root) {
	pi = (double)tot_pi/(double)num_procs;
	printf("%d trys pi = %.17g\n", trys, pi);
    }    

    MPI_Finalize();
    return 0;
} 

