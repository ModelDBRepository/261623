
//#include "neuron.h"
#include "layer.h"
#include <string>

int parse(Param* param)
{
	//read in parameters from parameters.txt
	//parameters to read in:Tgk, Tth, Th0, B, Te1, Te2, Ti1, Ti2 - both excitatory and inhibitory
	//Tau_Release, Tau_Replenish
	//Spontaneous Release Probability
	//Mean and Standard Deviation of synaptic depression
	//Connectivity distance, Small World probability
	//NRmax maximum
	/* set default values in returned structure */
	int i;
	string term;
	
	param->ExTgk   = 20; param->ExB	  = 20;	param->ExTth   = 10;	param->ExTe1   = 0.2;	param->ExTe2   = 10;	param->ExTi1	  = 0.2;	param->ExTi2   = 20;
	param->InTgk   = 10; param->InB	  = 10;	param->InTth   = 15;	param->InTe1   = 0.2;	param->InTe2   = 20;	param->InTi1	  = 0.2;	param->InTi2   = 20;	
	param->TauRelease = 1; param->TauReplenish = 5000; param->SpontRelease = 0.0001; param->Neighborhood = 5;
	param->NRmax = 15;

	cin >> term >> param->ExTgk;
	cin >> term >> param->ExB;
	cin >> term >> param->ExTth;
	cin >> term >> param->ExTe1;
	cin >> term >> param->ExTe2;
	cin >> term >> param->ExTi1;
	cin >> term >> param->ExTi2;
	
	cin >> term >> param->InTgk;
	cin >> term >> param->InB;
	cin >> term >> param->InTth;
	cin >> term >> param->InTe1;
	cin >> term >> param->InTe2;
	cin >> term >> param->InTi1;
	cin >> term >> param->InTi2;

	cin >> term >> param->TauRelease;
	cin >> term >> param->TauReplenish;
	cin >> term >> param->SpontRelease;
	cin >> term >> param->Neighborhood;
	cin >> term >> param->NRmax;


  return 0;
}

int main(){
	
	//omp_set_num_threads(20);
	//initialize random number generator
	srand (time(NULL));
	//srand (1);
	clock_t begin = clock();
#ifndef WIN32
	srand48(time(NULL));
	//srand48(1);
#endif
	Param* param;
	parse(param);
	//Specify the layer
	Layer pyramidLayer(100,100);
	//Add neurons
	pyramidLayer.addNeurons(param);
	//Add connections
	pyramidLayer.addSynapse(0, param);
	//Start simulation
	pyramidLayer.stimulateNeurons();	
	clock_t end = clock();
	double time_secs = double(end - begin) / CLOCKS_PER_SEC;
	//cout << "Time in secs:" << " " << time_secs << endl;
}