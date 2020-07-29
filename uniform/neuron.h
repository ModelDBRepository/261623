#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <time.h>
#include <omp.h>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

#ifdef WIN32
	#define drand48() ((double)rand()/(RAND_MAX));
#endif

//extern double factorial(double no);
//global variable, to scale the glutamate initial distribution
extern double distribution_max;
struct Param {
		double ExTgk,ExB,ExTth,ExTe1,ExTe2,ExTi1,ExTi2;
		double InTgk,InB,InTth,InTe1,InTe2,InTi1,InTi2;
		double TauRelease, TauReplenish, SpontRelease, Neighborhood;
		int NRmax;

};

struct Synapse{
	//indicate if the synapse is excitatory
	bool excitatory;
	//indices of presynaptic neuron
	int row, clm;
	//synaptic weight
	double weight;
	//vector of psps
	vector<double> psps;
	//glutamate release parameters:
	//NRmax: max releasable vesicles, NR: available releasable vesicles
	//PR: probability of release of releasable vesicle
	//tau: refill time constant
	double NRmax;
	double NR;
	vector <double> PR;
	double tau;
};

class Neuron{
	//parameters from the MacGregor Neuron Model. 
	//for details, refer: "A Model for Repetitive Firing in Neurons", MacGregor, R. J. and Oliver, R. M., Kybernetik, 1974.
	double Tmem, Th0;
	double mTmem, mTgk, mTh, mTCa;
	double dcGk, dcCa, Ge, Gi, GiGeRatio;
	double Vr, Ek, Ei, Ee;
	double c, CM, CS, Ca;	
	int P;
	vector<double> Ge_table, Gi_table, firing_Prob, glutamate_Prob;
public:
	double Tgk, Tth, Gk, B;
	double Te1, Te2, Ti1, Ti2; 
	double TauRelease, TauReplenish, SpontRelease;
	double E, Th, glutamate, gaba, glutamate_pr;
	int row, clm, S;
	bool excitatory;
	Neuron();
	int Inhibitory();
	int computePotential(vector<Synapse> &preSynapses, double SC, int time);
	int glutamateinitialize(vector<Synapse> &preSynapses,Neuron **NeuronLayer);
	int glutamaterelease(vector<Synapse> &preSynapses,Neuron **NeuronLayer, int iterations);
};
