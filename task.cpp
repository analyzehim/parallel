#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>

using namespace std;

int get_sign(int i, int j) //Get sign Hadamard matrix's element by indices
{
    int res = i & j;
    int count;
    for (count = 0; res; res >>= 1) 
        { 
            if (res & 1)
                count++;
        }
    if (count & 1 == 1)
        return -1;
    else
        return 1;
}

int print_HadamardMatrix(int N, int rank, int size) //print Hadamard matrix 2^n size
{

    int matr_size(pow(2,N));
    int *matr = new int[matr_size];

    int i_start = (matr_size / size) * rank;

    int i_end = (matr_size / size) * (rank + 1);

    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < matr_size; j++) 
        {

            if (get_sign(i,j) == 1) 
                matr[j] = 1;
            else 
                matr[j] = -1;
        }
    }
    return 0;
}

int main(int argc, char *argv[]) 
{
    int rank, size;

    MPI_Init (&argc, &argv);  /* starts MPI */
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);    /* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &size);    /* get number of processes */

    double t1 = MPI_Wtime();
    if(argc != 2) 
    {
        printf("Matrix size missing\n");
        return 1;
    }

    int N( atoi(argv[1]) );
    cout << "N: " << pow(2,N) << endl;

    print_HadamardMatrix(N, rank, size);

    double t2 = MPI_Wtime();
    cout << "TIME " << t2 - t1 << endl;
    MPI_Finalize();

    return 0;
}