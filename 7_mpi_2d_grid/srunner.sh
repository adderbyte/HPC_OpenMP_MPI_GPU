#!/bin/bash
#SBATCH --nodes 2


module load gcc mvapich2 scorep


time srun -n 12  ./gauss-seidel-mpi 60  60 4 3

