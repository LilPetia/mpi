#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 10 // Matrix size

void transpose(int *A, int *B, int n)
{
    //#pragma omp parallel for shared(A,B,n) // Parallel loop
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            B[i*n+j] = A[j*n+i]; 
        }
    }
}

int main()
{
    int *A = (int *)malloc(N * N * sizeof(int));
    int *B = (int *)malloc(N * N * sizeof(int));

    for (int i = 0; i < N*N; i++) {
        A[i] = rand() % 10;
    }
    // initialise the MPI 
    int ierr = MPI_Init(&argc, &argv);
    int procid, numprocs;

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if (ierr!= 0) {
        printf("MPI_Comm_rank failed\n");
        return 1;
    }
    for ()
    transpose(A, B, N);

    printf("Matrix A:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", A[i*N+j]);
        }
        printf("\n");
    }

    printf("Matrix B:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", B[i*N+j]);
        }
        printf("\n");
    }

    free(A);
    free(B);

    return 0;
}