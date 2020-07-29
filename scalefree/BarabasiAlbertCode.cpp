	int total_count = 0, temp_count = 0;
	//initial connection between the first 2 neurons
	i = 24; j = 24; m = 74; n = 74;
	count = 0;
	ConnectionMatrix[i][j].push_back(Synapse());
	ConnectionMatrix[i][j][count].row = m;
	ConnectionMatrix[i][j][count].clm = n;
	value = drand48();
	ConnectionMatrix[i][j][count].weight = value;
	ConnectionMatrix[i][j][count].excitatory = true;
	ConnectionMatrix[i][j][count].NRmax = 30*value;
	ConnectionMatrix[i][j][count].tau	= 0;
	No_of_Connections[i][j] += 1;
	ConnectionMatrix[m][n].push_back(Synapse());
	ConnectionMatrix[m][n][count].row = i;
	ConnectionMatrix[m][n][count].clm = j;
	value = drand48();
	ConnectionMatrix[m][n][count].weight = value;
	ConnectionMatrix[m][n][count].excitatory = true;
	ConnectionMatrix[m][n][count].NRmax = 30*value;
	ConnectionMatrix[m][n][count].tau	= 0;
	No_of_Connections[m][n] += 1;
	total_count += 2;
	//for(i = 0; i < rows; i++){
	//	for(j = 0; j < clms; j++){
	//		cout << No_of_Connections[i][j] << " " ;
	//	}
	//	cout << endl;
	//}
	//exit(0);
	//#pragma omp parallel for private(i,j,k,p,value,count,m,n)
	for(i = 0; i < rows; i++){
		for(j = 0; j < clms; j++){
			if(!((i == 24 && j == 24) || (i == 74 && j == 74))){
				count = 0;
				for(m = 0; m < rows; m++){
					for(n = 0; n < clms; n++){
						//inhibitory neurons do not see each other
						//no self feedback
						value = drand48();
						temp_count = total_count;
						if( (value < ((double)2*(double)No_of_Connections[m][n]/(double)temp_count)) && (!(i == m && j == n)) && (!(!NeuronLayer[i][j].excitatory && !NeuronLayer[m][n].excitatory)) ){
							ConnectionMatrix[i][j].push_back(Synapse());
							ConnectionMatrix[i][j][count].row = m;
							ConnectionMatrix[i][j][count].clm = n;
							value = drand48();
							ConnectionMatrix[i][j][count].weight = value;
							if(NeuronLayer[m][n].excitatory){
								ConnectionMatrix[i][j][count].excitatory = true;
								ConnectionMatrix[i][j][count].NRmax = 30*value;
								ConnectionMatrix[i][j][count].tau	= 0;
							}
							else{
								ConnectionMatrix[i][j][count].excitatory = false;
								ConnectionMatrix[i][j][count].NRmax = 0;
								ConnectionMatrix[i][j][count].tau	= 0;
								ConnectionMatrix[i][j][count].NR    = 0;
							}
							count = count+1;
							temp_count = ConnectionMatrix[m][n].size();
							ConnectionMatrix[m][n].push_back(Synapse());
							ConnectionMatrix[m][n][temp_count].row = i;
							ConnectionMatrix[m][n][temp_count].clm = j;
							value = drand48();
							ConnectionMatrix[m][n][temp_count].weight = value;
							if(NeuronLayer[i][j].excitatory){
								ConnectionMatrix[m][n][temp_count].excitatory = true;
								ConnectionMatrix[m][n][temp_count].NRmax = 30*value;
								ConnectionMatrix[m][n][temp_count].tau	= 0;
							}
							else{
								ConnectionMatrix[m][n][temp_count].excitatory = false;
								ConnectionMatrix[m][n][temp_count].NRmax = 0;
								ConnectionMatrix[m][n][temp_count].tau	= 0;
								ConnectionMatrix[m][n][temp_count].NR    = 0;
							}
							total_count += 2;
							No_of_Connections[i][j] += 1;
							No_of_Connections[m][n] += 1;
						}					
					}
				}
  			//cout << count << " " ;
			}
		}
	}
	cout << total_count << endl;
	for(i = 0; i < rows; i++){
		for(j = 0; j < clms; j++){
			cout << No_of_Connections[i][j] << " " ;
		}
		cout << endl;
	}

exit(0);	
