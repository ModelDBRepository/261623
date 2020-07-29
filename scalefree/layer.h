#include "neuron.h"


class Layer{
	int rows, clms;
	double **availableGlutamate, **releasableGlutamate;
	//NeuronLayer is a matrix of neurons
	Neuron **NeuronLayer;
	//ConnectionMatrix is a matrix of vectors; each vector contains the incoming synapses to the neuron at that location
	vector<Synapse> **ConnectionMatrix;
public:
	Layer(int row, int clm);
	int addNeurons(Param* param);
	int addSynapse(int no_wts, Param* param);
	int stimulateNeurons();
	int recordPotential(int iterations);
};
