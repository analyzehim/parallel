module add impi/4.1.0 slurm intel

mpicxx task.cpp

mpirun -n 1 task.cpp 10

sbatch -N 1 -n 5 -p gputest impi ./a.out 20