#!/bin/bash
#SBATCH --nodes 2
#SBATCH --ntasks 4

for N in  32 64 96 128 160 192 224 256 288 320
do
   ./gauss_correct  $N
done
