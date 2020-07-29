
#include "layer.h"

double distribution_max;

Layer::Layer(int row, int clm){
	rows = row;
	clms = clm;
	availableGlutamate = new double*[row]();
	for(int i = 0; i < row; i++)
		availableGlutamate[i] = new double[clm]();
	releasableGlutamate = new double*[row]();
	for(int i = 0; i < row; i++)
		releasableGlutamate[i] = new double[clm]();
}

int Layer::addNeurons(Param* param){
	int i, j;
	double mean, stddev;
	distribution_max = 0;
	mean = 50; stddev = 10.0;
	default_random_engine generator;
	normal_distribution<double> distribution(mean,stddev);
	//The initial plan was to pass the param pointer to the constructor. However, that would have involved initializing the neurons in a loop, as C++ does not have a method to initialize an array created with new without passing param to each instance in a loop.
	NeuronLayer = new Neuron*[rows]();
	for (i = 0; i < rows; i++){ 
		NeuronLayer[i] = new Neuron[clms]();
	}
	for(i = 0; i < rows; i++){
		for(j =0; j < clms; j++){
			availableGlutamate[i][j] = 0;
			releasableGlutamate[i][j] = 0;
			NeuronLayer[i][j].row = i;
			NeuronLayer[i][j].clm = j;		
			if((i%5 == 0) && (j%5 == 0)){
				NeuronLayer[i][j].excitatory = false;
				NeuronLayer[i][j].Inhibitory();
				NeuronLayer[i][j].Tgk = param->InTgk;
				NeuronLayer[i][j].Tth = param->InTth;
				NeuronLayer[i][j].B   = param->InB;
				NeuronLayer[i][j].Te1 = param->InTe1;
				NeuronLayer[i][j].Te2 = param->InTe2;
				NeuronLayer[i][j].Ti1 = param->InTi1;
				NeuronLayer[i][j].Ti2 = param->InTi2;
				//cout << i << " " << j << endl;
			}
			else{
				NeuronLayer[i][j].excitatory = true;
				NeuronLayer[i][j].Tgk = param->ExTgk;
				NeuronLayer[i][j].Tth = param->ExTth;
				NeuronLayer[i][j].B   = param->ExB;
				NeuronLayer[i][j].Te1 = param->ExTe1;
				NeuronLayer[i][j].Te2 = param->ExTe2;
				NeuronLayer[i][j].Ti1 = param->ExTi1;
				NeuronLayer[i][j].Ti2 = param->ExTi2;
				NeuronLayer[i][j].TauRelease	= param->TauRelease;
				NeuronLayer[i][j].TauReplenish	= param->TauReplenish;
				NeuronLayer[i][j].SpontRelease	= param->SpontRelease;
				NeuronLayer[i][j].glutamate_pr = distribution(generator);
				if(distribution_max < NeuronLayer[i][j].glutamate_pr)
					distribution_max = NeuronLayer[i][j].glutamate_pr;
			}
		}
	}
	return 0;
}

int Layer::addSynapse(int no_wts, Param* param){
	int i, j, k, m, n, count, connect, total_count = 0;
	double d = param->Neighborhood, value, p, in_Glut, wt;
	double **trackweight;
	trackweight = new double*[rows*clms]();
	for(int i = 0; i < rows*clms; i++)
		trackweight[i] = new double[rows*clms]();
	FILE *gmFile = fopen("glutamatemap.txt","w");
	FILE *igmFile = fopen("initialglutamatemap.txt","w");
	FILE *wtFile = fopen("weightmap.txt","w");
	ConnectionMatrix  = new vector<Synapse>*[rows]();
	for(i = 0; i < rows; i++){
		ConnectionMatrix[i] = new vector<Synapse>[clms]();
	}
//#pragma omp parallel for private(i,j,k,p,value,count,m,n)
	for(i = 0; i < rows; i++){
		for(j = 0; j < clms; j++){
			count = 0;
			in_Glut = 0;
			for(m = 0; m < rows; m++){
				for(n = 0; n < clms; n++){
					//initial connection criteria: distance between 2 neurons
					//no self feedback
					value = drand48();					
					if((value < exp(-((i - m)*(i - m))/9 -((j - n)*(j - n))/9)) && (sqrt((i - m)*(i - m) + (j - n)*(j - n)) < d) && (!(i == m && j == n)) && (!(!NeuronLayer[i][j].excitatory && !NeuronLayer[m][n].excitatory)) ){
					ConnectionMatrix[i][j].push_back(Synapse());
					ConnectionMatrix[i][j][count].row = m;
					ConnectionMatrix[i][j][count].clm = n;
					value = drand48();
					ConnectionMatrix[i][j][count].weight = value;
					trackweight[i*rows+j][m*rows+n] = value;
					if(NeuronLayer[m][n].excitatory){
						ConnectionMatrix[i][j][count].excitatory = true;
						ConnectionMatrix[i][j][count].NRmax = (param->NRmax)*value;
						ConnectionMatrix[i][j][count].tau	= 0;
					}
					else{
						ConnectionMatrix[i][j][count].excitatory = false;
						ConnectionMatrix[i][j][count].NRmax = 0;
						ConnectionMatrix[i][j][count].tau	= 0;
						ConnectionMatrix[i][j][count].NR    = 0;
					}
					in_Glut = in_Glut + ConnectionMatrix[i][j][count].NRmax;
					count = count+1;
					//fprintf(cFile,"%f ",value);
					}					
				}
			}
			//generating small world connectivity, following Watts-Strogatz model			
			for(k = 0; k < ConnectionMatrix[i][j].size(); k++){
				p = drand48();
				if( (p < 0.25) ){ 
					in_Glut = in_Glut - ConnectionMatrix[i][j][k].NRmax;
					trackweight[i*rows+j][ConnectionMatrix[i][j][k].row*rows+ConnectionMatrix[i][j][k].clm] = 0;
					ConnectionMatrix[i][j].erase(ConnectionMatrix[i][j].begin()+k);
					count = count - 1;					
					connect = 0;
					while(connect == 0){
						m = rows*drand48();
						n = clms*drand48();
						if((sqrt((i - m)*(i - m) + (j - n)*(j - n)) > d) && (!(i == m && j == n)) && (!(!NeuronLayer[i][j].excitatory && !NeuronLayer[m][n].excitatory)) ){
							connect = 1;
							ConnectionMatrix[i][j].push_back(Synapse());
							ConnectionMatrix[i][j][count].row = m;
							ConnectionMatrix[i][j][count].clm = n;
							value = drand48();
							ConnectionMatrix[i][j][count].weight = value;
							trackweight[i*rows+j][m*rows+n] = value;
							if(NeuronLayer[m][n].excitatory){
								ConnectionMatrix[i][j][count].excitatory = true;
								ConnectionMatrix[i][j][count].NRmax = (param->NRmax)*value;
								ConnectionMatrix[i][j][count].tau	= 0;
							}
							else{
								ConnectionMatrix[i][j][count].excitatory = false;
								ConnectionMatrix[i][j][count].NRmax = 0;
								ConnectionMatrix[i][j][count].tau	= 0;
								ConnectionMatrix[i][j][count].NR    = 0;
							}
							in_Glut = in_Glut + ConnectionMatrix[i][j][count].NRmax;
							count = count+1;
						}
					//fprintf(cFile,"%f ",0);
					}
				}
			}		
  			//cout << count << " " ;
			total_count = total_count + count;
			fprintf(gmFile,"%f ", in_Glut);
		}
		fprintf(gmFile,"\n");
	}
	//cout << total_count << endl;
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
			in_Glut = 0;
			for(count = 0; count < ConnectionMatrix[i][j].size(); count++){
				in_Glut = in_Glut + ConnectionMatrix[i][j][count].NR;
				if(ConnectionMatrix[i][j][count].excitatory)
					NeuronLayer[ConnectionMatrix[i][j][count].row][ConnectionMatrix[i][j][count].clm].glutamate += ConnectionMatrix[i][j][count].NR;
			}
			fprintf(igmFile,"%f ", in_Glut);
		}
		fprintf(igmFile,"\n");
	}
	for(i = 0; i < rows*clms; i++){
		for(j = 0; j < rows*clms; j++){
			fprintf(wtFile,"%f ", trackweight[i][j]);
		}
		fprintf(wtFile,"\n");
	}
	fclose(gmFile);
	fclose(igmFile);
	fclose(wtFile);
	return 0;
}

int Layer::stimulateNeurons(){
	int index, iterations;
	int i, j, espikes, ispikes, count;
    double avg_glutamate;
	FILE *pFile = fopen("spikes.txt","w");
	FILE *eFile = fopen("potential.txt","w");
	FILE *rrFile = fopen("releasedglutamate.txt","w");
	FILE *rgFile = fopen("releasableglutamate.txt","w");
	FILE *agFile = fopen("availableglutamate.txt","w");
	FILE *gbFile = fopen("gaba.txt","w");
	FILE *spFile = fopen("spikepatterns.txt","w");
	FILE *kFile = fopen("kconductance.txt","w");
	//variable index can be used to apply a stimulating current to random neuron(s) in each iteration
	for(iterations = 0; iterations < 60000; iterations++){
		index = rand()%(rows*clms);
		espikes = 0;
		ispikes = 0;
        avg_glutamate = 0;
		//loop to keep track of available glutamate vesicles at each neuron
		for(i = 0; i < rows; i++){
			for(j = 0; j < clms; j++){
				for(count = 0; count < ConnectionMatrix[i][j].size(); count++){
					if(ConnectionMatrix[i][j][count].excitatory){
						 availableGlutamate[ConnectionMatrix[i][j][count].row][ConnectionMatrix[i][j][count].clm] += ConnectionMatrix[i][j][count].NR;						
						}
					releasableGlutamate[i][j] += ConnectionMatrix[i][j][count].NR;//*ConnectionMatrix[i][j][count].weight;
					}
				}
		}
        //loop to record the average glutamate released and the number of excitatory and inhibitory neurons firing
		for(i = 0; i < rows; i++){
			for(j = 0; j < clms; j++){
              ///*if((iterations > 5000) && (iterations < 65000)){
					fprintf(eFile,"%f ", NeuronLayer[i][j].E);
					fprintf(rrFile,"%f ", NeuronLayer[i][j].glutamate);
					fprintf(rgFile,"%f ", releasableGlutamate[i][j]);
					fprintf(agFile,"%f ", availableGlutamate[i][j]);
					fprintf(gbFile,"%f ", NeuronLayer[i][j].gaba);
					fprintf(spFile,"%d ",NeuronLayer[i][j].S);
					fprintf(kFile,"%f ",NeuronLayer[i][j].Gk);
				//}*/
                avg_glutamate += NeuronLayer[i][j].glutamate;
				availableGlutamate[i][j] = 0;
				releasableGlutamate[i][j] = 0;
				if((i%5 == 0) && (j%5 == 0))
					ispikes += NeuronLayer[i][j].S;
				else
					espikes += NeuronLayer[i][j].S;
			}
		}
       ///*if((iterations > 5000) && (iterations < 65000)){
			fprintf(eFile,"\n");
			fprintf(rrFile,"\n");
			fprintf(rgFile,"\n");
			fprintf(agFile,"\n");
			fprintf(gbFile,"\n");
			fprintf(spFile,"\n");
			fprintf(kFile,"\n");
		//}*/
		//cout << espikes << " " << ispikes << " ";// << endl;
        //fprintf(gFile,"%f\n", avg_glutamate/(rows*clms));
		fprintf(pFile,"%d %d %d\n",iterations,espikes,ispikes);
#pragma omp parallel for private(i,j)
		for(i = 0; i < rows; i++){
			for(j = 0; j < clms; j++){				
				NeuronLayer[i][j].glutamaterelease(ConnectionMatrix[i][j],NeuronLayer,iterations);
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
	fclose(rrFile);
	fclose(rgFile);
	fclose(agFile);
	fclose(gbFile);
	fclose(spFile);
	fclose(kFile);
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
