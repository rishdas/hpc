/* C Example */
#include "mpi.h"
#include <math.h>
#include <stdio.h>

/*Reffered Code*/
/*https://www.dartmouth.edu/~rc/classes/intro_mpi/Numerical_integration_example.html*/
/*But the implemntation below is my own*/

float func(float x)
{
    return cos(x)*sin(x/2);
}

float integral(float start_point, int local_prob_size,
	       float width)
{
    float sample_width, start_point_n, result;

    result = 0.0;
    sample_width = width/2;
    for (int i = 0; i < local_prob_size; i++) {
        start_point_n = start_point + i*width;
        result += func(start_point_n + sample_width) * width;
    }
    return result;
}

int main(int argc,char *argv[])
{
    int        N = 6400000, local_prob_size;
    float      width, result, pi = acos(-1.0), lb = 0.0, ub = pi/2, result_l;
    float      start_point;
    float      local_width;

    int        rank, source, tag = 123, root = 0, num_procs;
    MPI_Status status;

    root = 0;
    tag = 123;


    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    width = (ub-lb)/N;
    local_prob_size = N/num_procs;
    local_width = (ub-lb)/num_procs;
    start_point = lb + rank*local_width;
    result_l = integral(start_point, local_prob_size, width);

    if(rank == root) {
        result = result_l;
        for (int i = 1; i < num_procs; i++) {
	    source = i;
	    MPI_Recv(&result_l, 1, MPI_REAL, source, tag,
		     MPI_COMM_WORLD, &status);
	    result += result_l;
        }
        printf("result = %f\n", result);
    } else {
        MPI_Send(&result_l, 1, MPI_REAL, root, tag,
		 MPI_COMM_WORLD);
    }
    MPI_Finalize();
}
