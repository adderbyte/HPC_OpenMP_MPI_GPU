#!/bin/bash
#SBATCH --nodes 1
for N in  32 64 96 128 160 192 224 256 288 320
do
   ./Siedel_execution  $N
done
