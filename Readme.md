##Summary
Hello! <Enter>
This is the works to the course "Supercomputer Modeling and Technology"

Author: Raev Evgeny, 611 group, CMC Faclulty

2016

## Task 1
The exercise was make the algowiki article with algorithm description.

My algorithm is generating Hadamar matrix.

Article available here: https://goo.gl/bgztFA

## Task 2
The exercise is the realise conjugate gradient method (also using steepest descent) for Puasson problem on uniform mesh.

Realisation via MPI & OpenMP.

##Useful commands on Lomonosov
module add impi/4.1.0 slurm intel<Enter>

mpicxx task.cpp<Enter>

mpirun -n 1 task.cpp 10<Enter>

sbatch -N 1 -n 5 -p gputest impi ./a.out 20<Enter>

scancel <job number>

###Have a nice day!