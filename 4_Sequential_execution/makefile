XX=gcc
LD=${CXX}
CXXFLAGS+=-Wall -Werror  -std=c++11 -lstdc++ 
LDFLAGS+=


Siedel_execution:Siedel_execution.o 
	$(LD) ${LDFLAGS} -o $@ $?
clean:
	rm -f hello *.o *~

