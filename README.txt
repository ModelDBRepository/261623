
The subfolder names uniform, smallworld, and scalefree refer to the network connectivity for the simulations contained in those folders.The code in the uniform folder has additional comments on the code, when compared to the other two. Each folder has their own READMEs as well. 

The name of the executable is macgregor. In each of the folders, you will also find the executables, scripts, and outputs for independent simulations with different values for the parameter of interest (neuron neighborhood for uniform, percentage of rewired connections in small world, and slope of line from scale free definitions for the scale free configuration).

Each of the folders contain: 1) Makefile 2) parameters.txt (used by the program to read in the various parameters) 3) main.cpp 4) neuron.cpp and neuron.h containing all of the code related to individual neurons 5) layer.cpp and layer.h containing all of the code related to the network layer. OpenMP is used to parallelize the for loops.

To create the executable, run the Makefile. To run the executable, do: ./macgregor < parameters.txt. The program writes out various parameters and outputs as can be seen in the code. The code for visualizing and analyzing the program output can be found in Matlab and Python folders of the root directory.
