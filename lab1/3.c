#include <mpi.h>
#include <stdio.h>

#define MSG_SIZE 1000
#define SAMPLES_SIZE 10000
int main(int argc, char *argv[]) {
  int rank, size;
  double t1, t2;
  char msg[MSG_SIZE];
  double samples_first_thread[SAMPLES_SIZE];
  double samples_second_thread[SAMPLES_SIZE];
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (size != 2) {
    printf("This program requires exactly 2 processes.\n");
    MPI_Finalize();
    return 1;
  }

  if (rank == 0) {
    for (int i = 0; i < SAMPLES_SIZE; i++) {
      t1 = MPI_Wtime();
      MPI_Send(msg, MSG_SIZE, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
      MPI_Recv(msg, MSG_SIZE, MPI_CHAR, 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      t2 = MPI_Wtime();
      samples_first_thread[i] = t2 - t1;
      printf("Communication time: %fs\n", t2 - t1);
    }
  } else {
    for (int i = 0; i < SAMPLES_SIZE; i++) {
      t1 = MPI_Wtime();
      MPI_Recv(msg, MSG_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Send(msg, MSG_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
      t2 = MPI_Wtime();
      samples_second_thread[i] = t2 - t1;
      printf("Communication time: %fs\n", t2 - t1);
    }
  }

  MPI_Finalize();
  FILE *fp1, *fp2;
  fp1 = fopen("output.csv", "w"); // create a file
  if (fp1 == NULL) {
    printf("Error while opening the file.\n");
    return 0;
  }
  for (size_t array_rowit = 0; array_rowit < SAMPLES_SIZE; array_rowit++) {
    fprintf(fp1, "%.10lf,", samples_first_thread[array_rowit]);
  }

  fclose(fp1);
  return 0;
}