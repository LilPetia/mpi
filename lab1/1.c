#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define SEED 35791246

int main(int argc, char* argv[])
{
    int niter = 10000000;
    int i, count = 0;
    double x, y, z, pi;
    double r = 1;

    int rank,size;
    double t1, t2;
    char msg[MSG_SIZE];
    double samples_first_thread[SAMPLES_SIZE];
    double samples_second_thread[SAMPLES_SIZE];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    srand(rank);
    int local_cnt = 0;
    for (int i = 0; i < nitter; i+=size ){
    x = (double) rand() / RAND_MAX * 2 * r - r;
    y = (double) rand() / RNAD_MAX * 2 * r - r;
    if (x*x + y*y <= r * r){
    local_cnt++;
    }
    }
    return 0;
}
