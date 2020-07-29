#include <iostream>
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

int factorial(int no)
{
  return (no == 1 || no == 0) ? 1 : factorial(no - 1) * no;
}

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
	double Tmem, Tgk, Tth, Th, Th0;
	double mTmem, mTgk, mTh, mTCa;
	double dcGk, dcCa, Ge, Gi, Gk, GiGeRatio;
	double Vr, Ek, Ei, Ee;
	double c, CM, CS, B, Ca;
	double Te1, Te2, Ti1, Ti2;  
	int P;
	vector<double> Ge_table, Gi_table, firing_Prob;
public:
	double E, glutamate;
	int row, clm, S;
	bool excitatory;
	Neuron();
	int Inhibitory();
	int computePotential(vector<Synapse> &preSynapses, double SC, int time);
	int glutamateinitialize(vector<Synapse> &preSynapses);
	int glutamaterelease(vector<Synapse> &preSynapses,Neuron **NeuronLayer);
};

Neuron::Neuron(){
	//initialize excitatory neurons
	Tmem  = 25, Tgk = 5, Tth = 25, Th = 0,Th0 = 1;
	mTmem = 1/Tmem, mTgk = 1/Tgk, mTh = 1/Tth;
	mTCa  = 1/270;
	dcGk  = exp(-mTgk), dcCa    = exp(-mTCa);
	Ge    = 0, Gi = 0, Gk    = 0.2, GiGeRatio = 3;
	Vr    = -70, Ek = -15, Ei = -15, Ee = 70, E = 0, Th = 0;
	c     = 0.75, CM = c, CS = 0.0, Ca = 0, B = 20; 
	//Te1	  = 0.2, Te2 = 0.3, Ti1 = 0.5, Ti2 = 2.0;  
	Te1	  = 0.3, Te2 = 1.5, Ti1 = 0.3, Ti2 = 1.5;  
	glutamate = 0;
	S	  = 0, P = 0;
	for(int i=0; i< 151; i++){
		Ge_table.push_back( (Te1 * Te2 / (Te1 - Te2)) * (exp(-i/Te1) - exp(-i/Te2)) );
		Gi_table.push_back( (Ti1 * Ti2 / (Ti1 - Ti2)) * (exp(-i/Ti1) - exp(-i/Ti2)) );
	}
	for(int i = 1; i <= 100; i++)
		firing_Prob.push_back( exp(-10000*(i/100 - 1)*(i/100 - 1)) );
}

int Neuron::Inhibitory(){
	//initialize inhibitory neurons
	Tmem  = 25, Tgk = 5, Tth = 25, Th = 0,Th0 = 1;
	mTmem = 1/Tmem, mTgk = 1/Tgk, mTh = 1/Tth;
	mTCa  = 1/270;
	dcGk  = exp(-mTgk), dcCa    = exp(-mTCa);
	Ge    = 0, Gi = 0, Gk    = 0.2;
	Vr    = -70, Ek = -10, Ei = -15, Ee = 70, E = 0, Th = 0;
	c     = 0.75, CM = c, CS = 0.0, Ca = 0, B = 20; 
	Te1	  = 0.3, Te2 = 1.5, Ti1 = 0.3, Ti2 = 1.5;  
	//Te1	  = 2, Te2 = 3, Ti1 = 5, Ti2 = 20; 
	S	  = 0;
	Ge_table.clear();
	Gi_table.clear();
	for(int i=0; i< 151; i++){
		Ge_table.push_back( (Te1 * Te2 / (Te1 - Te2)) * (exp(-i/Te1) - exp(-i/Te2)) );
		Gi_table.push_back( (Ti1 * Ti2 / (Ti1 - Ti2)) * (exp(-i/Ti1) - exp(-i/Ti2)) );
	}
	return 0;
}

int Neuron::computePotential(vector<Synapse> &preSynapses, double SC, int time){
	//compute potential of individual neurons
	unsigned int i,j;
	Ge = 0, Gi = 0, glutamate = 0;
	double mean = 0, stddev = 1.0, value;
	int prob_Index;
	default_random_engine generator;
	normal_distribution<double> distribution(mean,stddev);
	for(i = 0; i < preSynapses.size(); i++){
		//cout << preSynapses[i].psps.size() << " " << time << endl;
		j = preSynapses[i].psps.size();
		while(j >= 150){
			//preSynapses[i].psps.erase(preSynapses[i].psps.end()-1);
			preSynapses[i].psps.pop_back();
			//cout << "psp erase" << i << " " << j << endl;
			j = j-1;
		}
		if(preSynapses[i].excitatory){
			for(j = 0; j < preSynapses[i].psps.size(); j++){
				Ge += Ge_table[j]*preSynapses[i].psps[j];
				//to keep track of glutamate vesicles released before the neuron fires
				if(preSynapses[i].psps[j])
					glutamate += 1;
			}
		}
		else{
			for(j = 0; j < preSynapses[i].psps.size(); j++){
				Gi += Gi_table[j]*preSynapses[i].psps[j];
			}
		}
		/*while(j >= 100){
			preSynapses[i].psps.erase(preSynapses[i].psps.end()-1);
			//preSynapses[i].psps.pop_back();
			//cout << "psp erase" << i << " " << j << endl;
			j = j-1;
		}*/

	}
	//if((row%8 == 0) && (clm%8 == 0))
		//cout << "Gi-Ge-Gk" << " " << Gi << " " << Ge << " " << Gk << endl;	
	Gk = exp(-(double)1/Tgk)*Gk + B*S;
	Th = exp(-(double)1/Tth)*Th + (1 - exp(-(double)1/Tth))*(Th0+c*E);	
	//E  = exp(-(1+Gk+Ge+Gi)/5)*E + ((1 - exp(-(1+Gk+Ge+Gi)/5))*(SC+Ee*Ge+Ei*Gi+Ek*Gk+distribution(generator)))/(1+Gk+Ge+Gi);
	E  = exp(-((double)1+Gk+Ge+Gi)/(double)5)*E + (((double)1 - exp(-((double)1+Gk+Ge+Gi)/(double)5))*(SC+(double)7*Ge-(double)1*Gi-(double)1*Gk+distribution(generator)))/((double)1+Gk+Ge+Gi);
	//E  = exp(-(1+Gk+Ge+GiGeRatio)/25)*E + (1 - exp(-(1+Gk+Ge+GiGeRatio)/25))*(SC+ distribution(generator) + (Ee*Ge+Gi*Ei+Gk*Ek)/(1+Gk+Ge+GiGeRatio));//((Ge * layer->Ee + Gi * layer->Ei + Gk * layer->Ek) / Gtot + (noise + neuron->SC))
	if(E < -5)
		E = -5;
	if(E > 10)
		E = 10;
	
	if(S == 1){
		if(E < Th)
			S = 0;
	}
	if(S == 0){
            if((E > Th) && (P < 0)){
                S = 1;
                P = 3;
            }
		/*prob_Index = ceil(99*E/Th);
		//cout << prob_Index << endl;
		if(prob_Index < 0)
			prob_Index = 0;
		if(prob_Index > 99)
			prob_Index = 99;
		value = drand48();
		if(value < firing_Prob[prob_Index] && P < 0){
			S = 1;
			P = 4;
			//if(excitatory)
				//cout << glutamate << endl;
			//cout << "E" << " " << E << " " << "Th" << " " << Th << " " << S << endl;
		}*/

	}
	P = P-1;
    /*if(E >= Th && S == 0)
         S = 1;
     else
         S = 0;*/
	//if((row == 8) && (clm == 8))
		//cout << "E" << " " << E << " " << "Th" << " " << Th << " " << S << endl;

	


	return 0;
}

int Neuron::glutamateinitialize(vector<Synapse> &preSynapses){
	//initialize the number of releasable glutamate vesicles
	//the mean and stddev variables determine the initial glutamate distribution parameters
	int i, j;
	double p1, pr, nr;
	double mean, stddev, max;
	mean = 3.0; stddev = 5.0; max = 1.0;
	default_random_engine generator;
	normal_distribution<double> distribution(mean,stddev);
	
	for(i = 0; i < preSynapses.size();i++){
		if(preSynapses[i].excitatory){
			p1 = 1 - exp(-(double)1/(double)5);//(preSynapses[i].tau));
			nr = preSynapses[i].NRmax;
			pr  = 0;
			for(j = 0; j < nr; j++){
					pr += factorial(nr)*(pow(p1,j))*(pow((1 - p1),(nr-j)))/(factorial(j) * factorial(nr - j));
					//cout << pr << endl;
					preSynapses[i].PR.push_back(pr);	
			}
			preSynapses[i].NR = (distribution(generator)/max);
			if(preSynapses[i].NR  >= preSynapses[i].NRmax)
				preSynapses[i].NR  = preSynapses[i].NRmax - 1;
			if(preSynapses[i].NR  < 0)
				preSynapses[i].NR  = 0;
			//cout << "Glutamate Initialize" << " " << preSynapses[i].NR << endl;
		}
		else
			continue;
	}
	//exit(0);
	return 0;
}

int Neuron::glutamaterelease(vector<Synapse> &preSynapses,Neuron **NeuronLayer){
	int i;
	double p1, value;

	for(i = 0; i < preSynapses.size();i++){
		preSynapses[i].psps.insert(preSynapses[i].psps.begin(),0);
		if(preSynapses[i].excitatory){
			p1 = exp(-(double)1/(double)5000);//(preSynapses[i].tau);
			preSynapses[i].NR = preSynapses[i].NR*p1 + preSynapses[i].NRmax*(1 - p1);
			//cout << "NR" << " " << preSynapses[i].NR << " " << abs((int)ceil(preSynapses[i].NR)) << endl;

			if(preSynapses[i].NR  >= preSynapses[i].NRmax)
				preSynapses[i].NR  = preSynapses[i].NRmax - 1;
			if(preSynapses[i].NR  < 0)
				preSynapses[i].NR  = 0;
			//activity dependent glutamate release
			if(NeuronLayer[preSynapses[i].row][preSynapses[i].clm].S){
				value = drand48();
				if( value <= preSynapses[i].PR[(int)ceil(preSynapses[i].NR)] ){
					//glutamate release
					//cout << "psp add" << " " << value << " " << preSynapses[i].PR[(int)floor(preSynapses[i].NR)] << endl;
					value = preSynapses[i].weight;
					//preSynapses[i].psps.insert(preSynapses[i].psps.begin(),value);
					preSynapses[i].psps[0] = value;
					preSynapses[i].NR = preSynapses[i].NR - 1;
					if(preSynapses[i].NR  < 0)
						preSynapses[i].NR  = 0;
				}
			}
			//spontaneous release
			else{
				value = drand48();
				if(value < 0.001){
					value = preSynapses[i].weight;
					//preSynapses[i].psps.insert(preSynapses[i].psps.begin(),value);
					preSynapses[i].psps[0] = value;
					preSynapses[i].NR = preSynapses[i].NR - 1;
					if(preSynapses[i].NR  < 0)
						preSynapses[i].NR  = 0;
				}
			}
		}
		else{
			if(NeuronLayer[preSynapses[i].row][preSynapses[i].clm].S){
				value = drand48();
				if( value < 1){
					//GABA release
					//cout << "psp add" << " " << value << endl;
					value = preSynapses[i].weight;//drand48();
					//preSynapses[i].psps.insert(preSynapses[i].psps.begin(),value);
					preSynapses[i].psps[0] = value;
				}

			}
		}
		
	}
	
	//exit(0);
	return 0;
}

class Layer{
	int rows, clms;
	double **trackGlutamate;
	//NeuronLayer is a matrix of neurons
	Neuron **NeuronLayer;
	//ConnectionMatrix is a matrix of vectors; each vector contains the incoming synapses to the neuron at that location
	vector<Synapse> **ConnectionMatrix;
public:
	Layer(int row, int clm);
	int addNeurons();
	int addSynapse(int no_wts);
	int stimulateNeurons();
	int recordPotential(int iterations);
};

Layer::Layer(int row, int clm){
	rows = row;
	clms = clm;
	trackGlutamate = new double*[row]();
	for(int i = 0; i < row; i++)
		trackGlutamate[i] = new double[clm]();
}

int Layer::addNeurons(){
	int i, j;
	NeuronLayer = new Neuron*[rows]();
	for (i = 0; i < rows; i++){ 
		NeuronLayer[i] = new Neuron[clms]();
	}
	for(i = 0; i < rows; i++){
		for(j =0; j < clms; j++){
			NeuronLayer[i][j].row = i;
			NeuronLayer[i][j].clm = j;
			if((i%5 == 0) && (j%5 == 0)){
				NeuronLayer[i][j].excitatory = false;
				NeuronLayer[i][j].Inhibitory();
				//cout << i << " " << j << endl;
			}
			else
				NeuronLayer[i][j].excitatory = true;
		}
	}
	return 0;
}

int Layer::addSynapse(int no_wts){
	int i, j, m, n, count;
	double d = 20, value, p;
	ConnectionMatrix  = new vector<Synapse>*[rows]();
	for(i = 0; i < rows; i++){
		ConnectionMatrix[i] = new vector<Synapse>[clms]();
	}
#pragma omp parallel for private(i,j,count,m,n)
	for(i = 0; i < rows; i++){
		for(j = 0; j < clms; j++){
			count = 0;
			for(m = 0; m < rows; m++){
				for(n = 0; n < clms; n++){
                                
				//formation of a random network
				//no self feedback
                                        p = drand48();
					if((p < 0.5) && (sqrt((i - m)*(i - m) + (j - n)*(j - n)) < d) && (!(i == m && j == n)) && (!(!NeuronLayer[i][j].excitatory && !NeuronLayer[m][n].excitatory)) ){
						ConnectionMatrix[i][j].push_back(Synapse());
						ConnectionMatrix[i][j][count].row = m;
						ConnectionMatrix[i][j][count].clm = n;
						if(NeuronLayer[m][n].excitatory){
							ConnectionMatrix[i][j][count].excitatory = true;
							ConnectionMatrix[i][j][count].NRmax = 6;
							ConnectionMatrix[i][j][count].tau	= 1000;
						}
						else{
							ConnectionMatrix[i][j][count].excitatory = false;
							ConnectionMatrix[i][j][count].NRmax = 0;
							ConnectionMatrix[i][j][count].tau	= 0;
							ConnectionMatrix[i][j][count].NR    = 0;
						}
						value = drand48();
						ConnectionMatrix[i][j][count].weight = value;
						count = count+1;

					}
				}
			}
			//cout << count << " " ;
		}
	}
	//exit(0);
//#pragma omp parallel for private(i,j)
	for(i = 0; i < rows; i++){
		for(j = 0; j < clms; j++){
			NeuronLayer[i][j].glutamateinitialize(ConnectionMatrix[i][j]);
		}
	}
	for(i = 0; i < rows; i++){
		for(j = 0; j < clms; j++){
			for(count = 0; count < ConnectionMatrix[i][j].size(); count++){
				if(ConnectionMatrix[i][j][count].excitatory)
					NeuronLayer[ConnectionMatrix[i][j][count].row][ConnectionMatrix[i][j][count].clm].glutamate += ConnectionMatrix[i][j][count].NR;
			}
		}
	}
	return 0;
}

int Layer::stimulateNeurons(){
	int index, iterations;
	int i, j, espikes, ispikes;
        double avg_glutamate;
	FILE *pFile = fopen("spikes.txt","w");
	FILE *eFile = fopen("potential.txt","w");
	FILE *gFile = fopen("glutamate.txt","w");
	//variable index can be used to apply a stimulating current to random neuron(s) in each iteration
	for(iterations = 0; iterations < 10000; iterations++){
		index = rand()%(rows*clms);
		espikes = 0;
		ispikes = 0;
                avg_glutamate = 0;
                //loop to record the average glutamate released and the number of excitatory and inhibitory neurons firing
		for(i = 0; i < rows; i++){
			for(j = 0; j < clms; j++){
                            if((iterations > 5000) && (iterations < 6000))
                                fprintf(eFile,"%f ", NeuronLayer[i][j].E);
				//fprintf(gFile,"%f ", NeuronLayer[i][j].glutamate);
                                avg_glutamate += NeuronLayer[i][j].glutamate;
				if((i%5 == 0) && (j%5 == 0))
					ispikes += NeuronLayer[i][j].S;
				else
					espikes += NeuronLayer[i][j].S;
			}
		}
                if((iterations > 5000) && (iterations < 6000))
                    fprintf(eFile,"\n");
		//fprintf(gFile,"\n");
		//cout << iterations << " " << espikes << " " << ispikes << endl;
                fprintf(gFile,"%f\n", avg_glutamate/(rows*clms));
		fprintf(pFile,"%d %d %d\n",iterations,espikes,ispikes);
#pragma omp parallel for private(i,j)
		for(i = 0; i < rows; i++){
			for(j = 0; j < clms; j++){				
				NeuronLayer[i][j].glutamaterelease(ConnectionMatrix[i][j],NeuronLayer);
			}
		}
#pragma omp barrier
#pragma omp parallel for private(i,j)
		for(i = 0; i < rows; i++){
			for(j = 0; j < clms; j++){				
				//if(i*clms+j == index)
					NeuronLayer[i][j].computePotential(ConnectionMatrix[i][j],0, iterations);
				//else
					//NeuronLayer[i][j].computePotential(ConnectionMatrix[i][j],0, iterations);
			}
		}
#pragma omp barrier
		recordPotential(iterations);
	}
	fclose(eFile);
	fclose(gFile);
	fclose(pFile);
	return 0;
}

int Layer::recordPotential(int iterations){
	int i,j, center_i, center_j;
	double distance, potential, sigma = 3;
	center_i = rows/2;
	center_j = clms/2;
	//potential is assumed to be recorded from the center of the layer. 
	//neuron membrane potentials are multiplied by a Gaussian centered at the center.
	potential = 0.0;
	for(i = 0; i < rows; i++){
		for(j = 0; j < clms; j++){
			distance = ((double)1/sqrt((double)2*M_PI*sigma*sigma))*exp(-((center_i - i)*(center_i - i) + (center_j - j)*(center_j - j))/(double)2*sigma*sigma);
			potential += distance*NeuronLayer[i][j].E;
		}
	}
	cout << potential << " " << iterations << endl;
	return 0;
}

int main(){
	
	//omp_set_num_threads(20);

	//initialize random number generator
	srand (time(NULL));
	clock_t begin = clock();
#ifndef WIN32
	srand48(time(NULL));
#endif
	Layer pyramidLayer(100,100);
	pyramidLayer.addNeurons();
	pyramidLayer.addSynapse(0);
	pyramidLayer.stimulateNeurons();	
	clock_t end = clock();
	double time_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "Time in secs:" << " " << time_secs << endl;
}