#include<stdio.h>
#include<assert.h>
#include<stdlib.h>
#include<mpi.h>

typedef struct complex{
  double real;
  double imag;
} Complex;

int cal_pixel(Complex c){
  int count, max_iter;
  Complex z;
  double temp, lengthsq;
  
  max_iter = 256;
  z.real = 0;
  z.imag = 0;
  count = 0;

  do{
    temp = z.real * z.real - z.imag * z.imag + c.real;
    z.imag = 2 * z.real * z.imag + c.imag;
    z.real = temp;
    lengthsq = z.real * z.real + z.imag * z.imag;
    count ++;
  }
  while ((lengthsq < 4.0) && (count < max_iter));
  return(count);
}

#define MASTERPE 0
int main(int argc, char **argv){
  FILE *file;
  int i, j;
  Complex c;
  int tmp;
  double *data_l, *data_l_tmp;
  int nx, ny;
  int mystrt, myend;
  int nrows_l;
  int nprocs, mype;

  MPI_Status status;

  /* regular MPI initialization stuff */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  if (argc != 3){
    int err = 0;
    printf("argc %d\n", argc);
    if (mype == MASTERPE){
      printf("usage: mandelbrot nx ny");
      MPI_Abort(MPI_COMM_WORLD,err );
    }
  }

  /* get command line args */
  nx = atoi(argv[1]);
  ny = atoi(argv[2]);

  /* assume divides equally */
  nrows_l = nx/nprocs; 

  /* create buffer for local work only */
  data_l = (double *) malloc(nrows_l * ny * sizeof(double));
  data_l_tmp = data_l;

  /* calculate each processor's region of work */
  mystrt = mype*nrows_l;
  myend  = mystrt + nrows_l - 1;
  
  /* calc each procs coordinates and call local mandelbrot set function */
  for (i = mystrt; i <= myend; ++i){
    c.real = i/((double) nx) * 4. - 2. ; 
    for (j = 0; j < ny; ++j){
      c.imag = j/((double) ny) * 4. - 2. ; 
      tmp = cal_pixel(c);
      *data_l++ = (double) tmp;
    }
  }
  data_l = data_l_tmp;

  if (mype == MASTERPE){
    file = fopen("mandelbrot.bin_0000", "w");
    printf("nrows_l, ny  %d %d\n", nrows_l, ny);
    fwrite(data_l, nrows_l*ny, sizeof(double), file);
    fclose(file);

    for (i = 1; i < nprocs; ++i){
      MPI_Recv(data_l, nrows_l * ny, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      printf("received message from proc %d\n", i);
      file = fopen("mandelbrot.bin_0000", "a");
      fwrite(data_l, nrows_l*ny, sizeof(double), file);
      fclose(file);
    }
  }

  else{
    MPI_Send(data_l, nrows_l * ny, MPI_DOUBLE, MASTERPE, 0, MPI_COMM_WORLD);
  }  

  MPI_Finalize();
}

