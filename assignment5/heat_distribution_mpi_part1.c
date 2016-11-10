#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mpi.h"

#define COLS 12
#define ROWS 12
#define TEMP 50.
#define DEBUG 1
#define EPS 1e-3
#define I_FIX 5
#define J_FIX 5
#define WORKERS 4
#define ROOT 0
#define NONE -1
#define BEGIN_TAG 1
#define LTAG 2
#define RTAG 3
#define DONE 4
#define NOTDONE 5

double max_abs(double** m1, double** m2,
	       int start, int end){
    double max_val = DBL_MIN;
    for (int i = start; i <= end; i++)
        for (int j = 0; j < COLS; j++){
            if (fabs(m1[i][j] - m2[i][j]) > max_val) {
                max_val = fabs(m1[i][j] - m2[i][j]);
            }
        }
    return max_val;
}

void print_matrix(double** matrix, int start, int end){
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++)
            printf("%f ", matrix[i][j]);
        printf("\n");
    }
}

void copy_matrix(double** dest, double** source,
		 int start, int end) {
    for (int i = start; i <= end; i++)
        for (int j = 0; j < COLS; j++)
            dest[i][j] = source[i][j];
}

double** alloc_matrix(){
    double** matrix;
    matrix = (double**) malloc(ROWS * sizeof(double *));
    matrix[0] = (double*) malloc(ROWS * COLS * sizeof(double));
    for (int i = 1; i < ROWS; i++)
        matrix[i] = matrix[0] + i*COLS;
    return matrix;
}

void compute_new_values(double** old_matrix, double** new_matrix,
			int start, int end){
    for (int i = start; i <= end; i++)
        for (int j= 1; j < COLS-1; j++)
            new_matrix[i][j] =
                    0.25 * (old_matrix[i-1][j] + old_matrix[i+1][j] +
                            old_matrix[i][j-1] + old_matrix[i][j+1]);
    new_matrix[I_FIX][J_FIX] = TEMP;
}

void init_matrix(double** matrix, int rank){
    for (int i = 0; i < ROWS; i++)
        for (int j = 0; j < COLS; j++) {
            matrix[i][j] = 0.;
        }
    if (rank == ROOT) {
	matrix[I_FIX][J_FIX] = TEMP;
    }
}

int main(int argc, char *argv[]) {

    int        rank;
    int        no_of_workers;
    double     **a_old = alloc_matrix(); //allocate memory for the matrices
    double     **a_new = alloc_matrix();
    int        local_prob_size;
    int        offset;
    int        left, right;
    int        dest;
    int        source;
    int        msgtype;
    int        start, end;
    int        extra;
    MPI_Status status;
    int        rows;
    int        dummy = NOTDONE;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &no_of_workers);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    init_matrix(a_old, rank); //initialize the matrices
    init_matrix(a_new, rank);
    no_of_workers--; /*Root is not a worker rather a coordinator*/
    if (rank == ROOT) {
	while (1) {
	    local_prob_size = ROWS/no_of_workers;
	    extra = ROWS%no_of_workers;
	    offset = 0;
	    /*send all sub problems to workers*/
	    for (int i = 1; i<=no_of_workers; i++) {
		rows = (i<=extra) ? local_prob_size+1:local_prob_size;
		if (i == 1) {
		    left = NONE;
		} else {
		    left = i - 1;
		}
		if (i == no_of_workers) {
		    right = NONE;
		} else {
		    right = i + 1;
		}
		dest = i;
		/* printf("About to Sent to task %d: rows= %d offset= %d ", dest, */
		/* 	   local_prob_size, offset); */
		/* printf("left= %d right= %d\n", left, right); */

		MPI_Send(&offset, 1, MPI_INT, dest, BEGIN_TAG, MPI_COMM_WORLD);
		MPI_Send(&rows, 1, MPI_INT,
			 dest, BEGIN_TAG, MPI_COMM_WORLD);
		MPI_Send(&left, 1, MPI_INT,
			 dest, BEGIN_TAG, MPI_COMM_WORLD);
		MPI_Send(&right, 1, MPI_INT, dest, BEGIN_TAG, MPI_COMM_WORLD);
		MPI_Send(&a_new[offset][0], local_prob_size*COLS, MPI_DOUBLE,
			 dest, BEGIN_TAG, MPI_COMM_WORLD);
		printf("Sent to task %d: rows= %d offset= %d ", dest,
		       rows, offset);
		printf("left= %d right= %d\n", left, right);
		offset = offset + rows;
	    }
	    copy_matrix(a_old, a_new, 0, ROWS-1);
	    /* Now wait for results from all worker tasks */
	    for (int i=1; i<=no_of_workers-1; i++) {
		source = i;
		msgtype = DONE;
		MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, 
			 &status);
		MPI_Recv(&no_of_workers, 1, MPI_INT,
			 source, msgtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&a_new[offset][0], local_prob_size*COLS, MPI_DOUBLE, source,
			 msgtype, MPI_COMM_WORLD, &status);
	    }
	    double max_diff = max_abs(a_old, a_new, 0, ROWS-1);
	    if (DEBUG)
		printf("Max diff is: %f\n", max_diff);

	    if (max_diff < EPS) {
		dummy = DONE;
		MPI_Send(&dummy, 1, MPI_INT, dest, DONE, MPI_COMM_WORLD);
		break;
	    }
	}
    }
    if (rank != ROOT) {
	while(1) {
	    source = ROOT;
	    msgtype = BEGIN_TAG;
	    MPI_Recv(&offset, 1, MPI_INT, source, msgtype,
		     MPI_COMM_WORLD, &status);
	    MPI_Recv(&rows, 1, MPI_INT, source,
		     msgtype, MPI_COMM_WORLD, &status);
	    MPI_Recv(&left, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
	    MPI_Recv(&right, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
	    MPI_Recv(&a_new[offset][0], local_prob_size*COLS, MPI_DOUBLE,
		     source, msgtype, MPI_COMM_WORLD, &status);
	
	    start = offset;
	    end = offset + rows - 1;
	    if (offset==0) {
		start=1;
	    }
	    if ((offset+rows)==ROWS) {
		end--;
	    }
	    printf("task=%d  start=%d  end=%d\n", rank, start, end);


	    if (left != NONE) {
		MPI_Send(&a_new[offset][0], COLS, MPI_DOUBLE, left,
			 RTAG, MPI_COMM_WORLD);
		source = left;
		msgtype = LTAG;
		MPI_Recv(&a_new[offset-1][0], COLS, MPI_DOUBLE, source,
			 msgtype, MPI_COMM_WORLD, &status);
	    }
	    if (right != NONE) {
		MPI_Send(&a_new[end][0], COLS, MPI_DOUBLE, right,
			 LTAG, MPI_COMM_WORLD);
		source = right;
		msgtype = RTAG;
		MPI_Recv(&a_new[offset+rows][0], COLS, MPI_DOUBLE,
			 source, msgtype, MPI_COMM_WORLD, &status);
	    }
	    compute_new_values(a_old, a_new, start, end);
	    if (DEBUG) {
		printf("a_old is:\n"); //output matrix to screen
		print_matrix(a_old, start, end);

		printf("a_new is:\n");
		print_matrix(a_new, start, end);
	    }
	    /* double max_diff = max_abs(a_old, a_new, start, end); */
	    /* if (DEBUG) */
	    /*     printf("Max diff is: %f\n", max_diff); */

	    /* if (max_diff < EPS) */
	    /*     break; */

	    copy_matrix(a_old, a_new, start, end); //assign values of a_new to a_old

	    if (DEBUG)
		printf("End of iteration\n\n");
	    
	    MPI_Send(&offset, 1, MPI_INT, ROOT, DONE, MPI_COMM_WORLD);
	    MPI_Send(&local_prob_size, 1, MPI_INT, ROOT, DONE, MPI_COMM_WORLD);
	    MPI_Send(&a_new[offset][0], local_prob_size*COLS,
		     MPI_DOUBLE, ROOT, DONE, MPI_COMM_WORLD);
	    MPI_Recv(&dummy, 1, MPI_INT, source, msgtype,
		     MPI_COMM_WORLD, &status);
	    if (dummy == DONE) {
		break;
	    }
	}
	MPI_Finalize();
    }
    MPI_Finalize();
    printf("\nThe final heat distribution matrix is:\n");
    print_matrix(a_new, 0, ROWS);	

    /* printf("\nThe final heat distribution matrix is:\n"); */
    /* print_matrix(a_new); */

    return 0;
}
