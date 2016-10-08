#include "mpi.h"
#include <stdio.h>

typedef struct {
    int max_iter;
    double t0;
    double tf;
    double xmax[12];
    double xmin;
} Pars;

int main (int argc, char *argv[])
{
    int myid, numprocs, left, right;
    Pars buffer, buffer2;
    MPI_Request request;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    right = (myid + 1) % numprocs;
    left = (myid - 1);
    if (left < 0)
	left = numprocs - 1;

    //initiatize the send buffer
    buffer.max_iter = myid;
    buffer.t0 = 3.14*myid;
    buffer.tf = 1.67*myid;
    buffer.xmin = 2.55*myid;
    for  (int i = 0; i< 12; i++) {
	buffer.xmax[i] = 2.7*myid;
    }

    //Modify this
    //send myid to the left
    //MPI_sendrecv(&buffer, 1, FIXME, left, 123,
    //   &buffer2, 1, FIXME, right, 123, MPI_COMM_WORLD, &status);

    printf(" Process %d recieved %d\n", myid, buffer2.max_iter);

    MPI_Finalize();
    return 0;
}
