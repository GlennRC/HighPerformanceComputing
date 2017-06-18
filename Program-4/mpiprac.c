#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define ROOT 0
#define SIZE 256

// Creates an array of random numbers. Each number has a value from 0 - 1
void create_rand_nums(float** nums, int num_elements) {
  
  *nums = (float *)malloc(sizeof(float) * num_elements);
  
  for (int i = 0; i < num_elements; i++) {
    (*nums)[i] = (rand() / (float)RAND_MAX);
  }
}

int main(int argc, char** argv) {
    
    int comm_size, comm_rank, rc;

    rc = MPI_Init(&argc, &argv);

    if(rc != MPI_SUCCESS)
        MPI_Abort(MPI_COMM_WORLD, rc);

    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
    
    float* rand_nums;
    
    if(comm_rank == ROOT) {
        create_rand_nums(&rand_nums, SIZE);
        /**
        for(int i=0; i<SIZE; i++)
            printf("[%d] %f ", comm_rank, rand_nums[i]);
        printf("\n");
        **/
    }
        
    
    int proc_size = SIZE/comm_size;
    float* sub_nums = (float*)malloc(sizeof(float)*proc_size);

    MPI_Scatter(rand_nums, proc_size, MPI_FLOAT, sub_nums,
                    proc_size, MPI_FLOAT, ROOT, MPI_COMM_WORLD);

    float sum = 0;
    for(int i=0; i<proc_size; i++)
        sum += sub_nums[i];
    
    printf("[%d] sum = %f\n", comm_rank, sum);

    float* sums = NULL;
    
    if(comm_rank == ROOT)
        sums = (float*)malloc(sizeof(float)*comm_size);

    MPI_Gather(&sum, 1, MPI_FLOAT, sums, 1,
                    MPI_FLOAT, ROOT, MPI_COMM_WORLD);

    if(comm_rank == ROOT) {
        for(int i=0; i<comm_size; i++)
            printf("[%d] %f ", comm_rank, sums[i]);
        printf("\n");
    }
    
    MPI_Finalize();

    return 0;
}
