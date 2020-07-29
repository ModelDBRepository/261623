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

double factorial(double no)
{
  return ((no == 1 || no == 0) ? 1 : factorial(no - 1) * no);
}

double distribution_max;
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
	vector<double> Ge_table, Gi_table, firing_Prob, glutamate_Prob;
public:
	double E, glutamate, glutamate_pr;
	int row, clm, S;
	bool excitatory;
	Neuron();
	int Inhibitory();
	int computePotential(vector<Synapse> &preSynapses, double SC, int time);
	int glutamateinitialize(vector<Synapse> &preSynapses,Neuron **NeuronLayer);
	int glutamaterelease(vector<Synapse> &preSynapses,Neuron **NeuronLayer);
};

Neuron::Neuron(){
	//initialize excitatory neurons
	Tmem  = 25, Tgk = 20, Tth = 10, Th = 0,Th0 = 1;
	mTmem = 1/Tmem, mTgk = 1/Tgk, mTh = 1/Tth;
	mTCa  = 1/270;
	dcGk  = exp(-mTgk), dcCa    = exp(-mTCa);
	Ge    = 0, Gi = 0, Gk    = 1, GiGeRatio = 3;
	Vr    = -70, Ek = -15, Ei = -15, Ee = 70, E = 0, Th = 0;
	c     = 0.75, CM = c, CS = 0.0, Ca = 0, B = 20; 
	//Te1	  = 0.3, Te2 = 1.5, Ti1 = 0.3, Ti2 = 1.5;  
	Te1	  = 0.2, Te2 = 5, Ti1 = 0.2, Ti2 = 10;  
	glutamate = 0;
	S	  = 0, P = 0;
	for(int i=0; i< 151; i++){
		Ge_table.push_back( (Te1 * Te2 / (Te1 - Te2)) * (exp(-i/Te1) - exp(-i/Te2)) );
		Gi_table.push_back( (Ti1 * Ti2 / (Ti1 - Ti2)) * (exp(-i/Ti1) - exp(-i/Ti2)) );
	}
	for(int i = 1; i <= 100; i++)
		firing_Prob.push_back( exp(-1000*(i/100 - 1)*(i/100 - 1)) );
        for(int i = 1; i <= 3; i++)
		glutamate_Prob.push_back( exp(-1000*(i/4 - 1)*(i/4 - 1)) );
}

int Neuron::Inhibitory(){
	//initialize inhibitory neurons
	Tmem  = 25, Tgk = 10, Tth = 15, Th = 0,Th0 = 1;
	mTmem = 1/Tmem, mTgk = 1/Tgk, mTh = 1/Tth;
	mTCa  = 1/270;
	dcGk  = exp(-mTgk), dcCa    = exp(-mTCa);
	Ge    = 0, Gi = 0, Gk    = 1;
	Vr    = -70, Ek = -10, Ei = -15, Ee = 70, E = 0, Th = 0;
	c     = 0.75, CM = c, CS = 0.0, Ca = 0, B = 10; 
	//Te1	  = 0.3, Te2 = 1.5, Ti1 = 0.3, Ti2 = 1.5;  
	Te1	  = 0.2, Te2 = 10, Ti1 = 0.2, Ti2 = 10;  
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
	double mean = 0, stddev = 1.0, value, dist;
	int prob_Index;
	
	default_random_engine generator;
	normal_distribution<double> distribution(mean,stddev);
	for(i = 0; i < preSynapses.size(); i++){
		//cout << preSynapses[i].psps.size() << " " << time << endl;
		dist = sqrt((row - preSynapses[i].row)*(row - preSynapses[i].row) + (clm - preSynapses[i].clm)*(clm - preSynapses[i].clm));
		if(dist > 25)
			dist = 25;
		j = preSynapses[i].psps.size();
		while(j >= 150){
			preSynapses[i].psps.erase(preSynapses[i].psps.end()-1);
			//preSynapses[i].psps.pop_back();
			//cout << "psp erase" << i << " " << j << endl;
			j = j-1;
		}
		if(preSynapses[i].excitatory){
			for(j = round(dist); j < preSynapses[i].psps.size(); j++){
				Ge += Ge_table[j-round(dist)]*preSynapses[i].psps[j];
				//to keep track of glutamate vesicles released before the neuron fires
				if(preSynapses[i].psps[j])
					glutamate += 1;
			}
		}
		else{
			for(j = round(dist); j < preSynapses[i].psps.size(); j++){
				Gi += Gi_table[j - round(dist)]*preSynapses[i].psps[j];
			}
		}
		/*while(j >= 100){
			preSynapses[i].psps.erase(preSynapses[i].psps.end()-1);
			//preSynapses[i].psps.pop_back();
			//cout << "psp erase" << i << " " << j << endl;
			j = j-1;
		}*/

	}
	//if((row == 8) && (clm == 8))
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
            /*if((E > Th) && (P < 0)){
                S = 1;
                P = 3;
            }*/
		prob_Index = ceil(99*E/Th);
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
		}

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

int Neuron::glutamateinitialize(vector<Synapse> &preSynapses,Neuron **NeuronLayer){
	//initialize the number of releasable glutamate vesicles
	//the mean and stddev variables determine the initial glutamate distribution parameters
	double i, j, k;
	double p1, pr, nr;
	if(distribution_max == 0)
		distribution_max = 1;
	//double norm = 1;
	//for(i = 0; i < preSynapses.size();i++){
	//	norm = norm + preSynapses[i].weight * preSynapses[i].weight;
	//}
	//norm = pow(norm, 0.5);
	
	for(i = 0; i < preSynapses.size();i++){
		/*previous code*/
		if(preSynapses[i].excitatory){
			/*p1 = 1 - exp(-(double)1/(double)2);//(preSynapses[i].tau));
			//cout << p1 << endl;
			nr = preSynapses[i].NRmax;
			//pr  = 0;
			//preSynapses[i].PR.push_back(pr);	
			for(j = 1; j <= nr; j++){
				pr = 0;
				for( k = 1; k < j; k++){
					pr = pr + factorial(nr)*(pow(p1,k))*(pow((1 - p1),(nr-k)))/(factorial(k) * factorial(nr - k));
				}
				//cout << pr << endl;
				preSynapses[i].PR.push_back(pr);	
			}*/
			preSynapses[i].NR = preSynapses[i].NRmax*NeuronLayer[preSynapses[i].row][preSynapses[i].clm].glutamate_pr/distribution_max;
			if(preSynapses[i].NR  >= preSynapses[i].NRmax)
				preSynapses[i].NR  = preSynapses[i].NRmax;
			if(preSynapses[i].NR  < 0)
				preSynapses[i].NR  = 0;
			//cout << preSynapses[i].NR << endl;
		}
		else
			continue;
	}
	//exit(0);
	return 0;
}

int Neuron::glutamaterelease(vector<Synapse> &preSynapses,Neuron **NeuronLayer){
	int i,j,P_index;
	double p1, p2, value, Glu_Rel;
	p2 = exp(-(double)1/(double)1);

	for(i = 0; i < preSynapses.size();i++){
		preSynapses[i].psps.insert(preSynapses[i].psps.begin(),0);
		if(preSynapses[i].excitatory){
			p1 = exp(-(double)1/(double)5000);//(preSynapses[i].tau);
			preSynapses[i].NR = preSynapses[i].NR*p1 + preSynapses[i].NRmax*(1 - p1);
			preSynapses[i].tau += 1;
			if(preSynapses[i].NR  >= preSynapses[i].NRmax)
				preSynapses[i].NR  = preSynapses[i].NRmax;
			if(preSynapses[i].NR  < 0)
				preSynapses[i].NR  = 0;
			//activity dependent glutamate release
			if(NeuronLayer[preSynapses[i].row][preSynapses[i].clm].S){
                value = drand48();
				if( value <= 1){//PR[preSynapses[i].tau]){//[(int)round(preSynapses[i].NR)] ){
					//glutamate release					
					Glu_Rel = round(preSynapses[i].NR*p2);
					preSynapses[i].NR = preSynapses[i].NR - Glu_Rel;
					for(j = 0; j < Glu_Rel; j++){
						value = preSynapses[i].weight;
						preSynapses[i].psps[0] = value;
						if(j != Glu_Rel)
							preSynapses[i].psps.insert(preSynapses[i].psps.begin(),0);
					}
					////cout << "psp add" << " " << value << " " << preSynapses[i].PR[(int)floor(preSynapses[i].NR)] << endl;
					//value = preSynapses[i].weight;
					////preSynapses[i].psps.insert(preSynapses[i].psps.begin(),value);
					//preSynapses[i].psps[0] = value;
					//preSynapses[i].NR = preSynapses[i].NR - 1;
					preSynapses[i].tau = 0;
					if(preSynapses[i].NR  < 0)
						preSynapses[i].NR  = 0;
				}
                            //}
			}
			//spontaneous release
			else{
				value = drand48();
				if(value < 0.0001){
                    value = drand48();
                    if( value <= 1){//PR[preSynapses[i].tau]){//[(int)round(preSynapses[i].NR)] ){
					
					Glu_Rel = round(preSynapses[i].NR*p2);
					preSynapses[i].NR = preSynapses[i].NR - Glu_Rel;
					for(j = 0; j < Glu_Rel; j++){
						value = preSynapses[i].weight;
						preSynapses[i].psps[0] = value;
						if(j != Glu_Rel)
							preSynapses[i].psps.insert(preSynapses[i].psps.begin(),0);
					}

					//value = preSynapses[i].weight;
					////preSynapses[i].psps.insert(preSynapses[i].psps.begin(),value);
					//preSynapses[i].psps[0] = value;
					//preSynapses[i].NR = preSynapses[i].NR - 1;
					preSynapses[i].tau = 0;
					if(preSynapses[i].NR  < 0)
						preSynapses[i].NR  = 0;
                    }
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
			else{
				value = drand48();
				if( value < 0.001){
					//GABA spontaneous release
					//cout << "psp add" << " " << value << endl;
					//value = preSynapses[i].weight;//drand48();
					//preSynapses[i].psps.insert(preSynapses[i].psps.begin(),value);
					//preSynapses[i].psps[0] = value;
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
	double mean, stddev;
	distribution_max = 0;
	mean = 50; stddev = 10.0;
	default_random_engine generator;
	normal_distribution<double> distribution(mean,stddev);

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
				NeuronLayer[i][j].glutamate_pr = distribution(generator);
				if(distribution_max < NeuronLayer[i][j].glutamate_pr)
					distribution_max = NeuronLayer[i][j].glutamate_pr;
		}
	}
	return 0;
}

int Layer::addSynapse(int no_wts){
	int i, j, k, m, n, count, value_new, connect;
	double d = 5, value, p;
	
	std::default_random_engine generator;
	std::exponential_distribution<double> distribution(1000);
	ConnectionMatrix  = new vector<Synapse>*[rows]();
	for(i = 0; i < rows; i++){
		ConnectionMatrix[i] = new vector<Synapse>[clms]();
	}
	int **No_of_Connections;
	No_of_Connections = new int*[rows]();
	for(i = 0; i < rows; i++){
		No_of_Connections[i] = new int[clms]();
	}
	//379851//new: 364884
	//initially create a network where every neuron is connected to 20-40 neighbors.
	int total_count = 0, temp_count = 0;
	for(i = 0; i < rows; i++){
		for(j = 0; j < clms; j++){
			count = 0;
			value_new = floor(1000*distribution(generator));
			//if(value_new > 99)
				//value_new = 99;
			for(m = 0; m < rows; m++){
				for(n = 0; n < clms; n++){
					//initial connection criteria: distance between 2 neurons
					//no self feedback
					//&& (sqrt((i - m)*(i - m) + (j - n)*(j - n)) < d)
					//(value < exp(-((double)(i - m)*(i - m))/(double)connections[value_new] -((double)(j - n)*(j - n))/(double)connections[value_new]))
					value = drand48();					
					if((sqrt((i - m)*(i - m) + (j - n)*(j - n)) < (double)(d+value_new)) && (!(i == m && j == n)) && (!(!NeuronLayer[i][j].excitatory && !NeuronLayer[m][n].excitatory)) ){
					ConnectionMatrix[i][j].push_back(Synapse());
					ConnectionMatrix[i][j][count].row = m;
					ConnectionMatrix[i][j][count].clm = n;
					value = drand48();
					ConnectionMatrix[i][j][count].weight = value;
					if(NeuronLayer[m][n].excitatory){
						ConnectionMatrix[i][j][count].excitatory = true;
						ConnectionMatrix[i][j][count].NRmax = 15*value;
						ConnectionMatrix[i][j][count].tau	= 0;
					}
					else{
						ConnectionMatrix[i][j][count].excitatory = false;
						ConnectionMatrix[i][j][count].NRmax = 0;
						ConnectionMatrix[i][j][count].tau	= 0;
						ConnectionMatrix[i][j][count].NR    = 0;
					}
					count = count+1;
					//fprintf(cFile,"%f ",value);
					}					
				}
			}
  			//cout << count << " " ;
			total_count += count;
			No_of_Connections[i][j] += count;
		}
	}
	//cout << total_count << endl;
	//exit(0);	
	////initial connection between the first 2 neurons
	//i = 24; j = 24; m = 74; n = 74;
	//count = 0;
	//ConnectionMatrix[i][j].push_back(Synapse());
	//ConnectionMatrix[i][j][count].row = m;
	//ConnectionMatrix[i][j][count].clm = n;
	//value = drand48();
	//ConnectionMatrix[i][j][count].weight = value;
	//ConnectionMatrix[i][j][count].excitatory = true;
	//ConnectionMatrix[i][j][count].NRmax = 30*value;
	//ConnectionMatrix[i][j][count].tau	= 0;
	//No_of_Connections[i][j] += 1;
	//ConnectionMatrix[m][n].push_back(Synapse());
	//ConnectionMatrix[m][n][count].row = i;
	//ConnectionMatrix[m][n][count].clm = j;
	//value = drand48();
	//ConnectionMatrix[m][n][count].weight = value;
	//ConnectionMatrix[m][n][count].excitatory = true;
	//ConnectionMatrix[m][n][count].NRmax = 30*value;
	//ConnectionMatrix[m][n][count].tau	= 0;
	//No_of_Connections[m][n] += 1;
	//total_count += 2;
	////for(i = 0; i < rows; i++){
	////	for(j = 0; j < clms; j++){
	////		cout << No_of_Connections[i][j] << " " ;
	////	}
	////	cout << endl;
	////}
	////exit(0);
	////#pragma omp parallel for private(i,j,k,p,value,count,m,n)
	//for(i = 0; i < rows; i++){
	//	for(j = 0; j < clms; j++){
	//		if(!((i == 24 && j == 24) || (i == 74 && j == 74))){
	//			count = 0;
	//			for(m = 0; m < rows; m++){
	//				for(n = 0; n < clms; n++){
	//					//inhibitory neurons do not see each other
	//					//no self feedback
	//					value = drand48();
	//					temp_count = total_count;
	//					if( (value < ((double)2*(double)No_of_Connections[m][n]/(double)temp_count)) && (!(i == m && j == n)) && (!(!NeuronLayer[i][j].excitatory && !NeuronLayer[m][n].excitatory)) ){
	//						ConnectionMatrix[i][j].push_back(Synapse());
	//						ConnectionMatrix[i][j][count].row = m;
	//						ConnectionMatrix[i][j][count].clm = n;
	//						value = drand48();
	//						ConnectionMatrix[i][j][count].weight = value;
	//						if(NeuronLayer[m][n].excitatory){
	//							ConnectionMatrix[i][j][count].excitatory = true;
	//							ConnectionMatrix[i][j][count].NRmax = 30*value;
	//							ConnectionMatrix[i][j][count].tau	= 0;
	//						}
	//						else{
	//							ConnectionMatrix[i][j][count].excitatory = false;
	//							ConnectionMatrix[i][j][count].NRmax = 0;
	//							ConnectionMatrix[i][j][count].tau	= 0;
	//							ConnectionMatrix[i][j][count].NR    = 0;
	//						}
	//						count = count+1;
	//						temp_count = ConnectionMatrix[m][n].size();
	//						ConnectionMatrix[m][n].push_back(Synapse());
	//						ConnectionMatrix[m][n][temp_count].row = i;
	//						ConnectionMatrix[m][n][temp_count].clm = j;
	//						value = drand48();
	//						ConnectionMatrix[m][n][temp_count].weight = value;
	//						if(NeuronLayer[i][j].excitatory){
	//							ConnectionMatrix[m][n][temp_count].excitatory = true;
	//							ConnectionMatrix[m][n][temp_count].NRmax = 30*value;
	//							ConnectionMatrix[m][n][temp_count].tau	= 0;
	//						}
	//						else{
	//							ConnectionMatrix[m][n][temp_count].excitatory = false;
	//							ConnectionMatrix[m][n][temp_count].NRmax = 0;
	//							ConnectionMatrix[m][n][temp_count].tau	= 0;
	//							ConnectionMatrix[m][n][temp_count].NR    = 0;
	//						}
	//						total_count += 2;
	//						No_of_Connections[i][j] += 1;
	//						No_of_Connections[m][n] += 1;
	//					}					
	//				}
	//			}
 // 			//cout << count << " " ;
	//		}
	//	}
	//}
	//cout << total_count << endl;
//	for(i = 0; i < rows; i++){
//		for(j = 0; j < clms; j++){
//			cout << No_of_Connections[i][j] << " " ;
//		}
//		cout << endl;
//	}
//
//exit(0);	
//#pragma omp parallel for private(i,j)
	for(i = 0; i < rows; i++){
		for(j = 0; j < clms; j++){
			NeuronLayer[i][j].glutamateinitialize(ConnectionMatrix[i][j], NeuronLayer);
		}
	}
	//exit(0);	
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
	FILE *spFile = fopen("spikepatterns.txt","w");
	//variable index can be used to apply a stimulating current to random neuron(s) in each iteration
	for(iterations = 0; iterations < 60000; iterations++){
		index = rand()%(rows*clms);
		espikes = 0;
		ispikes = 0;
                avg_glutamate = 0;
                //loop to record the average glutamate released and the number of excitatory and inhibitory neurons firing
		for(i = 0; i < rows; i++){
			for(j = 0; j < clms; j++){
                //if((iterations > 5000) && (iterations < 15000))
				//fprintf(eFile,"%f ", NeuronLayer[i][j].E);
				//fprintf(gFile,"%f ", NeuronLayer[i][j].glutamate);
				//fprintf(spFile,"%f ",NeuronLayer[i][j].S);
                avg_glutamate += NeuronLayer[i][j].glutamate;
				if((i%5 == 0) && (j%5 == 0))
					ispikes += NeuronLayer[i][j].S;
				else
					espikes += NeuronLayer[i][j].S;
			}
		}
        //if((iterations > 5000) && (iterations < 15000))
        //fprintf(eFile,"\n");
		//fprintf(gFile,"\n");
		//fprintf(spFile,"\n");
		//cout << espikes << " " << ispikes << " ";// << endl;
                //fprintf(gFile,"%f\n", avg_glutamate/(rows*clms));
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
	fclose(spFile);
	fclose(pFile);
	return 0;
}

int Layer::recordPotential(int iterations){
	int i,j, center_i, center_j;
	double distance, potential, sigma = 2;
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
	//srand (1);
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
	//cout << "Time in secs:" << " " << time_secs << endl;
}