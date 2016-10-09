#include "mpi.h"
#include <stdio.h>
#include <stddef.h>

typedef struct {
    int max_iter;
    double t0;
    double tf;
    double xmax[12];
    double xmin;
} Pars;

MPI_Datatype create_mpi_pars_type();
int main (int argc, char *argv[])
{
    int myid, numprocs, left, right;
    Pars buffer, buffer2;
    MPI_Request request;
    MPI_Status status;
    MPI_Datatype mpi_pars_type;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    right = (myid + 2) % numprocs;
    left = (myid - 2);
    if (left < 0)
	left = numprocs + left;

    //initiatize the send buffer
    buffer.max_iter = myid;
    buffer.t0 = 3.14*myid;
    buffer.tf = 1.67*myid;
    buffer.xmin = 2.55*myid;
    for  (int i = 0; i< 12; i++) {
	buffer.xmax[i] = 2.7*myid;
    }

    mpi_pars_type = create_mpi_pars_type();
    //Modify this
    //send myid to the left
    MPI_Sendrecv(&buffer, 1, mpi_pars_type, left, 123,
      &buffer2, 1, mpi_pars_type, right, 123, MPI_COMM_WORLD, &status);

    printf(" Process %d recieved %d\n", myid, buffer2.max_iter);

    MPI_Finalize();
    return 0;
}
MPI_Datatype create_mpi_pars_type() {
    int nitems = 5;
    MPI_Datatype double_arr_type;
    MPI_Datatype struct_Pars_type;
    /*Create the sub type for double xmax*/
    MPI_Type_contiguous(12, MPI_DOUBLE, &double_arr_type);
    MPI_Type_commit(&double_arr_type);

    /*Creating the structure*/
    MPI_Datatype types[nitems];
    MPI_Aint offsets[nitems];
    int blocklengths[nitems];

    types[0] = MPI_INT; offsets[0] = offsetof(Pars, max_iter);
    blocklengths[0] = 1;

    types[1] = MPI_DOUBLE; offsets[1] = offsetof(Pars, t0);
    blocklengths[1] = 1;

    types[2] = MPI_DOUBLE; offsets[2] = offsetof(Pars, tf);
    blocklengths[2] = 1;

    types[3] = double_arr_type; offsets[3] = offsetof(Pars, xmax);
    blocklengths[3] = 12;

    types[4] = MPI_DOUBLE; offsets[4] = offsetof(Pars, xmin);
    blocklengths[4] = 1;

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &struct_Pars_type);
    MPI_Type_commit(&struct_Pars_type);

    return struct_Pars_type;

}
