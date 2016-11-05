#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define COLS 10
#define ROWS 10
#define TEMP 50.
#define DEBUG 1
#define EPS 1e-6
#define I_FIX 5
#define J_FIX 5

double max_abs(double** m1, double** m2){
    double max_val = DBL_MIN;
    for (int i = 0; i < COLS; i++)
        for (int j = 0; j < COLS; j++){
            if (fabs(m1[i][j] - m2[i][j]) > max_val) {
                max_val = fabs(m1[i][j] - m2[i][j]);
            }
        }
    return max_val;
}

void print_matrix(double** matrix){
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++)
            printf("%f ", matrix[i][j]);
        printf("\n");
    }
}

void copy_matrix(double** dest, double** source) {
    for (int i = 0; i < ROWS; i++)
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

void compute_new_values(double** old_matrix, double** new_matrix){
    for (int i = 1; i < ROWS-1; i++)
        for (int j= 1; j < COLS-1; j++)
            new_matrix[i][j] =
                    0.25 * (old_matrix[i-1][j] + old_matrix[i+1][j] +
                            old_matrix[i][j-1] + old_matrix[i][j+1]);
    new_matrix[I_FIX][J_FIX] = TEMP;
}

void init_matrix(double** matrix){
    for (int i = 0; i < ROWS; i++)
        for (int j = 0; j < COLS; j++) {
            matrix[i][j] = 0.;
        }
    matrix[I_FIX][J_FIX] = TEMP;
}

int main(int argc, char *argv[]) {

    double **a_old = alloc_matrix(); //allocate memory for the matrices
    double **a_new = alloc_matrix();

    init_matrix(a_old); //initialize the matrices
    init_matrix(a_new);

    while (1) {

        if (DEBUG)
            printf("Performing a new iteration...\n");

        //compute new values and put them into a_new
        compute_new_values(a_old, a_new);

        if (DEBUG) {
            printf("a_old is:\n"); //output matrix to screen
            print_matrix(a_old);

            printf("a_new is:\n");
            print_matrix(a_new);
        }

        //calculate the maximum absolute differences among pairwise
        // differences of old and new matrix elements
        double max_diff = max_abs(a_old, a_new);

        if (DEBUG)
            printf("Max diff is: %f\n", max_diff);

        if (max_diff < EPS)
            break;

        copy_matrix(a_old, a_new); //assign values of a_new to a_old

        if (DEBUG)
            printf("End of iteration\n\n");
    }

    printf("\nThe final heat distribution matrix is:\n");
    print_matrix(a_new);

    return 0;
}
