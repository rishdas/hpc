#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define COLS        10000
#define ROWS        10000
#define STEPS       10000
#define BEGIN       1
#define LTAG        2
#define RTAG        3
#define NONE        0
#define DONE        4
#define NOTDONE     5
#define ROOT        0
#define EPS         1e-3
#define I_FIX       5000
#define J_FIX       5000
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
void compute_new_values(double** old_matrix, double** new_matrix,
			  int start, int end) {
    
    for (int i = start; i <= end; i++) {

	for (int j= 1; j < COLS-1; j++) {

	    new_matrix[i][j] =
		0.25 * (old_matrix[i-1][j] + old_matrix[i+1][j] +
			old_matrix[i][j-1] + old_matrix[i][j+1]);
	}
    }

    new_matrix[I_FIX][J_FIX] = TEMP;
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
void printdat(int r, int c, double **mat, char *fnam) {
    int i, j;
    FILE *fp;

    fp = fopen(fnam, "w");
    for (i = I_FIX; i < (I_FIX+50); i++) {
	for (j = J_FIX; j < (J_FIX+50); j++) {
	    fprintf(fp, "%f", mat[i][j]);
	    if (j != ((J_FIX+50)-1)) 
		fprintf(fp, " ");
	    else
		fprintf(fp, "\n");
	}
    }
    fclose(fp);
}
int main (int argc, char *argv[])
{
    double      **a_old = alloc_matrix();
    double      **a_new = alloc_matrix();
    int         no_of_procs;
    int	        rank;
    int         no_of_workers;
    int         row_per_worker;
    int         rows,offset,extra;
    int         dest, source;
    int         left,right;
    int         msgtype;
    int         start,end;
    int         i;
    int         task_done = 0;
    int         all_done = 0;
    int         ctr = 0;
    double      max_diff = 0.;
    double      g_max_diff = 0.;
    MPI_Status  status;
    MPI_Request request;
    MPI_Request srrequest;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&no_of_procs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    no_of_workers = no_of_procs-1;
    
    if (rank == ROOT) {

	printf ("Starting %d worker procs\n", no_of_workers);

	init_matrix(a_old, rank);
	init_matrix(a_new, rank);
//	print_matrix(a_new);
	printdat(COLS, ROWS, a_old, "initial.dat");

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
	/* while(1) { */
	/*     /\* printf("rank: %d maxdiff: %f g_max_diff: %f\n", *\/ */
	/*     /\* 	   rank, max_diff, g_max_diff); *\/ */
	/*     max_diff = 0; */
	/*     MPI_Iallreduce(&max_diff, &g_max_diff, 1, */
	/* 		   MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, &srrequest); */
	/*     ctr++; */
	/*     if (ctr>STEPS && g_max_diff < EPS) { */
	/* 	break; */
	/*     } */
	/*     MPI_Wait(&srrequest, &status); */
	/* } */

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
	//print_matrix(a_new);
	printdat(COLS, ROWS, a_new, "final.dat");
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

	while (1) {
	    if (right != NONE) {
		MPI_Isend(&a_old[offset+rows-1][0], COLS, MPI_DOUBLE, right,
			  LTAG, MPI_COMM_WORLD, &srrequest);
		source = right;
		msgtype = RTAG;
		MPI_Irecv(&a_old[offset+rows][0], COLS, MPI_DOUBLE,
			  source, msgtype, MPI_COMM_WORLD, &srrequest);
	    }

	    if (left != NONE) {
		MPI_Isend(&a_old[offset][0], COLS, MPI_DOUBLE, left,
			  RTAG, MPI_COMM_WORLD, &srrequest);
		source = left;
		msgtype = LTAG;
		MPI_Irecv(&a_old[offset-1][0], COLS, MPI_DOUBLE, source,
			  msgtype, MPI_COMM_WORLD, &srrequest);
	    }

	    MPI_Wait(&srrequest, &status);
	    compute_new_values(a_old, a_new, start, end);	    
	    max_diff = max_abs(a_old, a_new, start, end);
	    //g_max_diff = max_diff;
	    copy_matrix(a_old, a_new, start, end);	    
	    printf("rank: %d maxdiff: %f g_max_diff: %f\n",
		   rank, max_diff, g_max_diff);

	    ctr++;
	    if (ctr>STEPS) {// && g_max_diff < EPS) {
		break;
	    }
	    MPI_Iallreduce(&max_diff, &g_max_diff, 1,
			   MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, &srrequest);
	}
	MPI_Send(&offset, 1, MPI_INT, ROOT, DONE, MPI_COMM_WORLD);
	MPI_Send(&rows, 1, MPI_INT, ROOT, DONE, MPI_COMM_WORLD);
	MPI_Send(&a_new[offset][0], rows*COLS, MPI_DOUBLE, ROOT, DONE, 
		 MPI_COMM_WORLD);
	MPI_Finalize();
    }
}

