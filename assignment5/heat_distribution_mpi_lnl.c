#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define COLS        20
#define ROWS        20
#define STEPS       100
#define BEGIN       1
#define LTAG        2
#define RTAG        3
#define NONE        0
#define DONE        4
#define NOTDONE     5
#define ROOT        0
#define EPS         1e-3
#define I_FIX       10
#define J_FIX       10
#define TEMP        50.

double** alloc_matrix() {
    
    double** matrix;
    matrix = (double**) malloc(ROWS * sizeof(double *));
    matrix[0] = (double*) malloc(ROWS * COLS * sizeof(double));
    for (int i = 1; i < ROWS; i++)
        matrix[i] = matrix[0] + i*COLS;
    return matrix;
    
}
void init_matrix(double** matrix, int rank) {
    
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            matrix[i][j] = 0.;
        }
    }

    if (rank == ROOT) {
	matrix[I_FIX][J_FIX] = TEMP;
    }
}
void copy_matrix(double** dest, double** source,
		 int start, int end) {
    
    for (int i = start; i <= end; i++) {

	for (int j = 0; j < COLS; j++) {
            dest[i][j] = source[i][j];
	}
    }
}
double compute_new_values(double** old_matrix, double** new_matrix,
			  int start, int end) {
    
    double max_val = DBL_MIN;
    for (int i = start; i <= end; i++) {

	for (int j= 1; j < COLS-1; j++) {

	    new_matrix[i][j] =
		0.25 * (old_matrix[i-1][j] + old_matrix[i+1][j] +
			old_matrix[i][j-1] + old_matrix[i][j+1]);
	    if (fabs(new_matrix[i][j] - old_matrix[i][j]) > max_val) {
                max_val = fabs(new_matrix[i][j] - old_matrix[i][j]);
            }
	}
    }

    new_matrix[I_FIX][J_FIX] = TEMP;
    return max_val;
}
void print_matrix(double** matrix) {
    
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            printf("%f ", matrix[i][j]);
	}
        printf("\n");
    }
}
void print_array(double** matrix, int row) {
    
    for (int j = 0; j < COLS; j++) {
	printf("%f ", matrix[row][j]);
    }

    printf("\n");
}
double max_abs(double** m1, double** m2,
	       int start, int end) {
    
    double max_val = DBL_MIN;
    for (int i = start; i <= end; i++) {
        for (int j = 0; j < COLS; j++){
            if (fabs(m1[i][j] - m2[i][j]) > max_val) {
                max_val = fabs(m1[i][j] - m2[i][j]);
            }
        }
    }
    return max_val;
}
int main (int argc, char *argv[])
{
    double     **a_old = alloc_matrix();
    double     **a_new = alloc_matrix();
    int        no_of_procs;
    int	       rank;
    int        no_of_workers;
    int        row_per_worker;
    int        rows,offset,extra;
    int        dest, source;
    int        left,right;
    int        msgtype;
    int        start,end;
    int        i;
    int        task_done = 0;
    int        all_done = 0;
    int        ctr = 0;
    double     max_diff = 0.;
    MPI_Status status;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&no_of_procs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    no_of_workers = no_of_procs-1;
    
    if (rank == ROOT) {

	printf ("Starting %d worker procs\n", no_of_workers);

	init_matrix(a_old, rank);
	init_matrix(a_new, rank);
	print_matrix(a_new);

	offset = 0;
		
	row_per_worker = ROWS/no_of_workers;
	extra = ROWS%no_of_workers;
	
	for (i=1; i<=no_of_workers; i++) {

	    rows = (i <= extra) ? row_per_worker+1 : row_per_worker; 

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

	    MPI_Send(&offset, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
	    MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
	    MPI_Send(&left, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
	    MPI_Send(&right, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
	    MPI_Send(&a_new[offset][0], rows*COLS, MPI_DOUBLE, dest, BEGIN, 
		     MPI_COMM_WORLD);

	    printf("Sent to rank %d: rows= %d offset= %d ",dest,rows,offset);
	    printf("left= %d right= %d\n",left,right);
	    offset = offset + rows;
	}
	copy_matrix(a_old, a_new, 0, ROWS-1);

	for (i=1; i<=no_of_workers; i++) {
	    source = i;
	    msgtype = DONE;
	    MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, 
		     &status);
	    MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
	    MPI_Recv(&a_new[offset][0], rows*COLS, MPI_DOUBLE, source,
		     msgtype, MPI_COMM_WORLD, &status);
	    printf("Recieving results from : %d\n", source);
	}
	print_matrix(a_new);	
	MPI_Finalize();
    }
    if (rank != ROOT) {
	source = ROOT;
	msgtype = BEGIN;
	MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&left, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&right, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&a_new[offset][0], rows*COLS, MPI_DOUBLE, source, msgtype, 
		 MPI_COMM_WORLD, &status);

	start=offset;
	end=offset+rows-1;

	if (offset==0)  {
	    start=1;
	}
	
	if ((offset+rows)==ROWS) {
	    end--;
	}
	
	printf("task=%d  start=%d  end=%d\n",rank,start,end);
	printf("Task %d received work. Beginning time steps...\n",rank);

	copy_matrix(a_old, a_new, start, end);
	
	while (ctr < STEPS) {
	    if (left != NONE)
	    {
		MPI_Send(&a_new[offset][0], COLS, MPI_DOUBLE, left,
			 RTAG, MPI_COMM_WORLD);
		source = left;
		msgtype = LTAG;
		MPI_Recv(&a_new[offset-1][0], COLS, MPI_DOUBLE, source,
			 msgtype, MPI_COMM_WORLD, &status);
	    }
	    if (right != NONE)
	    {
		MPI_Send(&a_new[offset+rows-1][0], COLS, MPI_DOUBLE, right,
			 LTAG, MPI_COMM_WORLD);
		source = right;
		msgtype = RTAG;
		MPI_Recv(&a_new[offset+rows][0], COLS, MPI_DOUBLE,
			 source, msgtype, MPI_COMM_WORLD, &status);
	    }
	
	    max_diff = compute_new_values(a_old, a_new, start, end);	    
	    /* max_diff = max_abs(a_old, a_new, start, end); */
	    if (max_diff < EPS) {
		task_done = 1;	       
	    }
	    copy_matrix(a_old, a_new, start, end);
	    //printf("rank: %d ALL DONE: %d\n", rank, all_done);
	    ctr++;
	}
	MPI_Send(&offset, 1, MPI_INT, ROOT, DONE, MPI_COMM_WORLD);
	MPI_Send(&rows, 1, MPI_INT, ROOT, DONE, MPI_COMM_WORLD);
	MPI_Send(&a_new[offset][0], rows*COLS, MPI_DOUBLE, ROOT, DONE, 
		 MPI_COMM_WORLD);
	MPI_Finalize();
    }
}

