EXECS=gauss_correct
CC=mpicc
CFLAGS=-Wall -Werror -g  
LDFLAGS=-lm 

all: ${EXECS}

gauss_correct:gauss_correct.cc utils.cc
	mpicc  -O3 -std=c++11 -lstdc++ gauss_correct.cc utils.cc -o  gauss_correct -lm 

clean:
	rm ${EXECS} 
