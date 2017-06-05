.PHONY=all
all: gauss-seidel-mpi
CC=mpicc
CFLAGS=-Wall -Werror -g
LDFLAGS=-lm
gauss-seidel-mpi: gauss-seidel-mpi.c utils.c
	scorep  mpicc -O3 -Wall -D TEST_CONV -D PRINT_RESULTS gauss-seidel-mpi.c utils.c -o gauss-seidel-mpi -lm

.PHONY: clean

