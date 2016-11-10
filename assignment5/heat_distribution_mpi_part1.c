#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define NXPROB      20                 /* x dimension of problem grid */
#define NYPROB      20                 /* y dimension of problem grid */
#define COLS 20
#define ROWS 20
#define STEPS       100                /* number of time steps */
#define MAXWORKER   8                  /* maximum number of worker tasks */
#define MINWORKER   3                  /* minimum number of worker tasks */
#define BEGIN       1                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define NONE        0                  /* indicates no neighbor */
#define DONE        4                  /* message tag */
#define NOTDONE     5                  /* message tag */
#define ROOT        0                  /* taskid of first process */
#define EPS         1e-3
#define I_FIX 5
#define J_FIX 5
#define TEMP 50.

double** alloc_matrix(){
    double** matrix;
    matrix = (double**) malloc(ROWS * sizeof(double *));
    matrix[0] = (double*) malloc(ROWS * COLS * sizeof(double));
    for (int i = 1; i < ROWS; i++)
        matrix[i] = matrix[0] + i*COLS;
    return matrix;
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
void copy_matrix(double** dest, double** source,
		 int start, int end) {
    for (int i = start; i <= end; i++)
        for (int j = 0; j < COLS; j++)
            dest[i][j] = source[i][j];
}
double compute_new_values(double** old_matrix, double** new_matrix,
			int start, int end){
    double max_val = DBL_MIN;
    for (int i = start; i <= end; i++)
        for (int j= 1; j < COLS-1; j++) {
            new_matrix[i][j] =
                    0.25 * (old_matrix[i-1][j] + old_matrix[i+1][j] +
                            old_matrix[i][j-1] + old_matrix[i][j+1]);
	    if (fabs(new_matrix[i][j] - old_matrix[i][j]) > max_val) {
                max_val = fabs(new_matrix[i][j] - old_matrix[i][j]);
            }
	}
    new_matrix[I_FIX][J_FIX] = TEMP;
    return max_val;
}
void print_matrix(double** matrix) {
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++)
            printf("%f ", matrix[i][j]);
        printf("\n");
    }
}
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
int main (int argc, char *argv[])
{
    void inidat(), prtdat(), update();
    int	taskid,                     /* this task's unique id */
	numworkers,                 /* number of worker processes */
	numtasks,                   /* number of tasks */
	averow,rows,offset,extra,   /* for sending rows of data */
	dest, source,               /* to - from for message send-receive */
	left,right,        /* neighbor tasks */
	msgtype,                    /* for message types */
	rc,start,end,               /* misc */
	i,ix,iy,iz,it;              /* loop variables */
    MPI_Status status;
    double **a_old = alloc_matrix(); //allocate memory for the matrices
    double **a_new = alloc_matrix();
    int    dummy;


/* First, find out my taskid and how many tasks are running */
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    numworkers = numtasks-1;

    init_matrix(a_old, taskid); //initialize the matrices
    init_matrix(a_new, taskid);

    if (taskid == ROOT) {
	/************************* master code *******************************/
	/* Check if numworkers is within range - quit if not */
	if ((numworkers > MAXWORKER) || (numworkers < MINWORKER)) {
	    printf("ERROR: the number of tasks must be between %d and %d.\n",
		   MINWORKER+1,MAXWORKER+1);
	    printf("Quitting...\n");
	    MPI_Abort(MPI_COMM_WORLD, rc);
	    exit(1);
	}
	printf ("Starting mpi_heat2D with %d worker tasks.\n", numworkers);

	/* Initialize grid */
	printf("Grid size: X= %d  Y= %d  Time steps= %d\n",NXPROB,NYPROB,STEPS);
	init_matrix(a_old, taskid); //initialize the matrices
	init_matrix(a_new, taskid);
	print_matrix(a_new);

	/* Distribute work to workers.  Must first figure out how many rows to */
	/* send and what to do with extra rows.  */
	averow = NXPROB/numworkers;
	extra = NXPROB%numworkers;
	offset = 0;
	for (i=1; i<=numworkers; i++)
	{
	    rows = (i <= extra) ? averow+1 : averow; 
	    /* Tell each worker who its neighbors are, since they must exchange */
	    /* data with each other. */  
	    if (i == 1) 
		left = NONE;
	    else
		left = i - 1;
	    if (i == numworkers)
		right = NONE;
	    else
		right = i + 1;
	    /*  Now send startup information to each worker  */
	    dest = i;
	    MPI_Send(&offset, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
	    MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
	    MPI_Send(&left, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
	    MPI_Send(&right, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
	    MPI_Send(&a_new[offset][0], rows*NYPROB, MPI_DOUBLE, dest, BEGIN, 
		     MPI_COMM_WORLD);
	    printf("Sent to task %d: rows= %d offset= %d ",dest,rows,offset);
	    printf("left= %d right= %d\n",left,right);
	    offset = offset + rows;
	}
	copy_matrix(a_old, a_new, 0, ROWS-1);
	/* Now wait for results from all worker tasks */
	for (i=1; i<=numworkers; i++)
	{
	    source = i;
	    msgtype = DONE;
	    MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, 
		     &status);
	    MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
	    MPI_Recv(&a_new[offset][0], rows*NYPROB, MPI_DOUBLE, source,
		     msgtype, MPI_COMM_WORLD, &status);
	    printf("Recieving results from : %d\n", source);
	}

	/* Write final output, call X graph and finalize MPI */
	print_matrix(a_new);	
	MPI_Finalize();
    }   /* End of master code */



    /************************* workers code **********************************/
    if (taskid != ROOT) 
    {
	/* Receive my offset, rows, neighbors and grid partition from master */
	source = ROOT;
	msgtype = BEGIN;
	MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&left, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&right, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
	MPI_Recv(&a_new[offset][0], rows*NYPROB, MPI_DOUBLE, source, msgtype, 
		 MPI_COMM_WORLD, &status);

	/* Determine border elements.  Need to consider first and last columns. */
	/* Obviously, row 0 can't exchange with row 0-1.  Likewise, the last */
	/* row can't exchange with last+1.  */
	start=offset;
	end=offset+rows-1;
	if (offset==0) 
	    start=1;
	if ((offset+rows)==NXPROB) 
	    end--;
	printf("task=%d  start=%d  end=%d\n",taskid,start,end);

	/* Begin doing STEPS iterations.  Must communicate border rows with */
	/* neighbors.  If I have the first or last grid row, then I only need */
	/*  to  communicate with one neighbor  */
	printf("Task %d received work. Beginning time steps...\n",taskid);
	copy_matrix(a_old, a_new, start, end);
	int ctr = 1;
	double max_diff = 0.;
	printf("start %d\n", taskid);
	while (ctr < STEPS) {
	    ctr++;
	    printf("Task: %d: Iteration %d: max_diff: %f EPS: %f\n",
	    	   taskid, ctr, max_diff, EPS);

	    if (left != NONE)
	    {
		MPI_Send(&a_new[offset][0], NYPROB, MPI_DOUBLE, left,
			 RTAG, MPI_COMM_WORLD);
		source = left;
		msgtype = LTAG;
		MPI_Recv(&a_new[offset-1][0], NYPROB, MPI_DOUBLE, source,
			 msgtype, MPI_COMM_WORLD, &status);
	    }
	    if (right != NONE)
	    {
		MPI_Send(&a_new[offset+rows-1][0], NYPROB, MPI_DOUBLE, right,
			 LTAG, MPI_COMM_WORLD);
		source = right;
		msgtype = RTAG;
		MPI_Recv(&a_new[offset+rows][0], NYPROB, MPI_DOUBLE,
			 source, msgtype, MPI_COMM_WORLD, &status);
	    }
	    /* Now call update to update the value of grid points */
	    //update(start,end,NYPROB,&a_new[0][0],&a_old[0][0]);
	    max_diff = compute_new_values(a_old, a_new, start, end);	    
	    /* max_diff = max_abs(a_old, a_new, start, end); */
	    /* if (max_diff < EPS) */
	    /* 	break; */
	    copy_matrix(a_old, a_new, start, end);	
	}
	/* Finally, send my portion of final results back to master */
	MPI_Send(&offset, 1, MPI_INT, ROOT, DONE, MPI_COMM_WORLD);
	MPI_Send(&rows, 1, MPI_INT, ROOT, DONE, MPI_COMM_WORLD);
	MPI_Send(&a_new[offset][0], rows*NYPROB, MPI_DOUBLE, ROOT, DONE, 
		 MPI_COMM_WORLD);
	MPI_Finalize();
    }
}


/**************************************************************************
 *  subroutine update
 ****************************************************************************/
void update(int start, int end, int ny, double *u1, double *u2)
{
    int ix, iy;
    for (ix = start; ix <= end; ix++) 
	for (iy = 1; iy <= ny-2; iy++) 
	    *(u2+ix*ny+iy) = 0.25 * (*(u1+(ix+1)*ny+iy) + *(u1+(ix-1)*ny+iy) +
				     *(u1+ix*ny+iy+1) + *(u1+ix*ny+iy-1));
}

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, double *u) {
    int ix, iy;

    for (ix = 0; ix <= nx-1; ix++) 
	for (iy = 0; iy <= ny-1; iy++)
	    *(u+ix*ny+iy) = 0.;

    *(u+I_FIX*ny+J_FIX) = TEMP;
}

/**************************************************************************
 * subroutine prtdat
 **************************************************************************/
void prtdat(int nx, int ny, double *u1, char *fnam) {
    int ix, iy;
    FILE *fp;

    fp = fopen(fnam, "w");
    for (iy = ny-1; iy >= 0; iy--) {
	for (ix = 0; ix <= nx-1; ix++) {
	    fprintf(fp, "%f", *(u1+ix*ny+iy));
	    if (ix != nx-1) 
		fprintf(fp, " ");
	    else
		fprintf(fp, "\n");
	}
    }
    fclose(fp);
}
