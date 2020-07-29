
#include "neuron.h"


Neuron::Neuron(){
	//initialize excitatory neurons
	Tmem  = 25, Tgk = 20, Tth = 15, Th = 0,Th0 = 1;
	mTmem = 1/Tmem, mTgk = 1/Tgk, mTh = 1/Tth;
	mTCa  = 1/270;
	dcGk  = exp(-mTgk), dcCa    = exp(-mTCa);
	Ge    = 0, Gi = 0, Gk    = 1, GiGeRatio = 3;
	Vr    = -70, Ek = -15, Ei = -15, Ee = 70, E = 0, Th = 0;
	c     = 0.75, CM = c, CS = 0.0, Ca = 0, B = 20; 
	//Te1	  = 0.3, Te2 = 1.5, Ti1 = 0.3, Ti2 = 1.5;  
	Te1	  = 0.2, Te2 = 10, Ti1 = 0.2, Ti2 = 20;  
	TauRelease = 1; TauReplenish = 8000; SpontRelease = 0.0001;
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
	Te1	  = 0.2, Te2 = 20, Ti1 = 0.2, Ti2 = 20;  
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
	Ge = 0, Gi = 0, glutamate = 0, gaba = 0;
	double mean = 0, stddev = 1.0, value, dist;
	int prob_Index;
	
	default_random_engine generator;
	normal_distribution<double> distribution(mean,stddev);
	for(i = 0; i < preSynapses.size(); i++){
		dist = sqrt((row - preSynapses[i].row)*(row - preSynapses[i].row) + (clm - preSynapses[i].clm)*(clm - preSynapses[i].clm));
		if(dist > 25)
			dist = 25;
		j = preSynapses[i].psps.size();
		while(j >= 150){
			preSynapses[i].psps.erase(preSynapses[i].psps.end()-1);
			j = j-1;
		}
		if(preSynapses[i].excitatory){
			for(j = round(dist); j < preSynapses[i].psps.size(); j++){
				Ge += Ge_table[j-round(dist)]*preSynapses[i].psps[j];
				//to keep track of glutamate vesicles released before the neuron fires
				if(preSynapses[i].psps[j])
					glutamate += 1;//preSynapses[i].psps[j];
			}
		}
		else{
			for(j = round(dist); j < preSynapses[i].psps.size(); j++){
				Gi += Gi_table[j - round(dist)]*preSynapses[i].psps[j];
				if(preSynapses[i].psps[j])
					gaba += 1;//preSynapses[i].psps[j];
			}
		}

	}

	Gk = exp(-(double)1/Tgk)*Gk + B*S;
	Th = exp(-(double)1/Tth)*Th + (1 - exp(-(double)1/Tth))*(Th0+c*E);	
	E  = exp(-((double)1+Gk+Ge+Gi)/(double)5)*E + (((double)1 - exp(-((double)1+Gk+Ge+Gi)/(double)5))*(SC+(double)7*Ge-(double)1*Gi-(double)1*Gk+distribution(generator)))/((double)1+Gk+Ge+Gi);
	if(E < -5)
		E = -5;
	if(E > 10)
		E = 10;
	
	if(S == 1){
		if(E < Th)
			S = 0;
	}
	if(S == 0){
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
		}

	}
	P = P-1;
	return 0;
}

int Neuron::glutamateinitialize(vector<Synapse> &preSynapses,Neuron **NeuronLayer){
	//initialize the number of releasable glutamate vesicles
	//the mean and stddev variables determine the initial glutamate distribution parameters
	double i, j, k;
	double p1, pr, nr;
	if(distribution_max == 0)
		distribution_max = 1;
	
	for(i = 0; i < preSynapses.size();i++){
		if(preSynapses[i].excitatory){
			preSynapses[i].NR = preSynapses[i].NRmax*NeuronLayer[preSynapses[i].row][preSynapses[i].clm].glutamate_pr/distribution_max;
			if(preSynapses[i].NR  >= preSynapses[i].NRmax)
				preSynapses[i].NR  = preSynapses[i].NRmax;
			if(preSynapses[i].NR  < 0)
				preSynapses[i].NR  = 0;
		}
		else
			continue;
	}
	return 0;
}

int Neuron::glutamaterelease(vector<Synapse> &preSynapses,Neuron **NeuronLayer,int iterations){
	int i,j,P_index;
	double p1, p2, value, Glu_Rel;
	p2 = exp(-(double)1/TauRelease);
	//glutamate = 0;
	for(i = 0; i < preSynapses.size();i++){
		preSynapses[i].psps.insert(preSynapses[i].psps.begin(),0);
		if(preSynapses[i].excitatory){
			p1 = exp(-(double)1/TauReplenish);//(preSynapses[i].tau);
			preSynapses[i].NR = preSynapses[i].NR*p1 + preSynapses[i].NRmax*(1 - p1);
			//glutamate += preSynapses[i].NR/preSynapses[i].NRmax ;
			preSynapses[i].tau += 1;
			if(preSynapses[i].NR  >= preSynapses[i].NRmax)
				preSynapses[i].NR  = preSynapses[i].NRmax;
			if(preSynapses[i].NR  < 0)
				preSynapses[i].NR  = 0;
			//activity dependent glutamate release
			if(NeuronLayer[preSynapses[i].row][preSynapses[i].clm].S){
                value = drand48();
				if( value <= 1 ){
					//glutamate release					
					Glu_Rel = round(preSynapses[i].NR*pow(p2,4 - NeuronLayer[preSynapses[i].row][preSynapses[i].clm].P));
					preSynapses[i].NR = preSynapses[i].NR - Glu_Rel;
					//glutamate = glutamate + Glu_Rel;
					for(j = 0; j < Glu_Rel; j++){
						value = preSynapses[i].weight;
						preSynapses[i].psps[0] = value;
						if(j != Glu_Rel)
							preSynapses[i].psps.insert(preSynapses[i].psps.begin(),0);
					}
					preSynapses[i].tau = 0;
					if(preSynapses[i].NR  < 0)
						preSynapses[i].NR  = 0;
				}
                            //}
			}
			//spontaneous release
			else{
				value = drand48();
				if(value < SpontRelease){ 
                    value = drand48();
                    if( value <= 1 ){ 
					Glu_Rel = round(preSynapses[i].NR*p2);
					preSynapses[i].NR = preSynapses[i].NR - Glu_Rel;
					//glutamate = glutamate + Glu_Rel;
					for(j = 0; j < Glu_Rel; j++){
						value = preSynapses[i].weight;
						preSynapses[i].psps[0] = value;
						if(j != Glu_Rel)
							preSynapses[i].psps.insert(preSynapses[i].psps.begin(),0);
					}

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
				if( value < 1 ){
					//GABA release
					//cout << "psp add" << " " << value << endl;
					value = preSynapses[i].weight;//drand48();
					preSynapses[i].psps[0] = value;
				}

			}
			else{
				value = drand48();
				if( value < SpontRelease){
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
