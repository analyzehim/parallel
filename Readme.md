#Summary
Hello! <Enter>

#Useful commands on Lomonosov
module add impi/4.1.0 slurm intel<Enter>

mpicxx task.cpp<Enter>

mpirun -n 1 task.cpp 10<Enter>

sbatch -N 1 -n 5 -p gputest impi ./a.out 20<Enter>