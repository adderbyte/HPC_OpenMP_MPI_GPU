#!/bin/bash
#SBATCH --nodes 1
#SBATCH --reservation phpc2017
#SBATCH --account phpc2017
for N in 32 64 96 128 160 192 224 256 288 320 352 384 416 448 480 512 544 576 608 640 672 704 736 768 800 832 864 896 928 960 992 1024
do
   ./Siedel_execution  $N
done
