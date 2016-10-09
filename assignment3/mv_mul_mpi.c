#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

/* The idea for the code written below has been taken from
 * https://people.sc.fsu.edu/~jburkardt/c_src/mpi/matvec_mpi.c
 * But the implementation has been written by me.
 */
int main(int argc, char *argv[])
{
    double *mat, *mat_row, *vect, *result;
    double res;
    int dummy;
    int i, j;
    int mat_index;
    int N = 1000;
    int root = 0;
    int rank;
    int num_procs;
    int row_no;
    int active_procs;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (rank == root) {
	mat = (double*)malloc(N*N*sizeof(double));
	vect = (double*)malloc(N*sizeof(double));
	result = (double*)malloc(N*sizeof(double));
	mat_index = 0;
	for ( i = 1; i <= N; i++ ) {
	    for ( j = 1; j <= N; j++ ) {
		mat[mat_index] = 1.3*i*j;
		mat_index = mat_index + 1;
	    }
	}
	for ( i = 0; i < N; i++ ) {
	    vect[i] = 2.7*i;
	}
    } else {
	mat_row = (double*)malloc(N*sizeof(double));
	vect = (double*)malloc(N*sizeof(double));
    }
    MPI_Bcast(vect, N, MPI_DOUBLE, root, MPI_COMM_WORLD);

    if (rank == root) {
	row_no = 0;
	for (i = 1; i <= num_procs-1; i++) {
	    mat_index = row_no*N;
	    MPI_Send(mat+mat_index, N, MPI_DOUBLE, i, row_no, MPI_COMM_WORLD);
	    row_no++;
	}
	active_procs = num_procs-1;
	while (1) {
	    MPI_Recv(&res, 1, MPI_DOUBLE, MPI_ANY_SOURCE,
		     MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	    result[status.MPI_TAG] = res;
	    if (row_no < N) {
		row_no++;
		mat_index = row_no*N;
		MPI_Send(mat+mat_index, N, MPI_DOUBLE, status.MPI_SOURCE,
			 row_no, MPI_COMM_WORLD);
	    } else {
		active_procs--;
		dummy = 0;
		MPI_Send(&dummy, 1, MPI_INT, status.MPI_SOURCE,
			 N+1, MPI_COMM_WORLD);
		if(active_procs == 0) {
		    break;
		}
	    }
	}
	free(mat);
	free(vect);
    } else {
	while (1) {
	    MPI_Recv(mat_row, N, MPI_DOUBLE, root, MPI_ANY_TAG,
		     MPI_COMM_WORLD, &status);
	    if (status.MPI_TAG == N+1) {
		break;
	    }
	    res = 0.0;
	    for (i = 0; i < N; i++) {
		res = res + mat_row[i]*vect[i];
	    }
	    MPI_Send(&res, 1, MPI_DOUBLE, root, status.MPI_TAG,
		     MPI_COMM_WORLD);
	}
	free(vect);
	free(mat_row);
    }
    if (rank == root) {
	for (i = 0; i < N; i++)
	    printf ( "%g\n", result[i]);
	free(result);
    }
    MPI_Finalize();
    return 0;
}
