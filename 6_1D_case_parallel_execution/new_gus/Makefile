.PHONY:all
all: gaus
CC=mpicc
CFLAGS=-Wall -Werror -g
LDFLAGS=-lm

gaus: gaus.c utils.c
	mpicc -O3 -Wall gaus.c utils.c -o gaus -lm

.PHONY:clean	
