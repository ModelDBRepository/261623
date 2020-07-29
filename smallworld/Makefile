COMPILER = g++
CCFLAGS = -std=c++11 -g -fopenmp 

macgregor:main.o layer.o neuron.o
	${COMPILER} ${CCFLAGS} -o macgregor main.o layer.o neuron.o
neuron.o: neuron.cpp neuron.h
	${COMPILER} ${CCFLAGS} -c neuron.cpp  
layer.o: layer.cpp layer.h
	${COMPILER} ${CCFLAGS} -c layer.cpp  
main.o: main.cpp layer.h 
	${COMPILER} ${CCFLAGS} -c main.cpp  
clean: 
	rm -rf *.o macgregor

