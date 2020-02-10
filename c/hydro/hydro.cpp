#include "database.h"
#include "nauty.h"
#include "hydro.h"
#include "sampling.h"
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include "../defines.h"

#define HYDRO_CUT 1.08
namespace bd {

// ********************************************************************** //
// ********************* Class Functions ******************************** //
// ********************************************************************** //

HydroCluster::HydroCluster() {
	//constructor

	cNum = 0;
	
}

HydroCluster::~HydroCluster() {
	//deconstructor

	delete []clusters;
}

HCC::HCC(int num_clusters_, int N_, int maxT_) {
	//constructor

	num_clusters = num_clusters_;
	N = N_;
	maxT = maxT_;

	hc = new HydroCluster[num_clusters];

	for (int i = 0; i < num_clusters; i++) {
		(hc)[i].N = N_; (hc)[i].maxT = maxT_; 
		if (N_ == 6) {
			(hc)[i].cNum = 1;
		}
		(hc)[i].clusters = new double[N_*DIMENSION*maxT_];
		for (int j = 0; j < N_* DIMENSION * maxT_; j++) {
			(hc)[i].clusters[j] = 0;
		}
	}
}

HCC::~HCC() {
	//deconstructor

	delete []hc;
}

HCC* extractData(std::string& filename, int N, int maxT) {
	//read data from HD file and put into a HC class

	std::ifstream in_str(filename);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		return NULL;
	}

	//read first line, N = number of particles
	int num_clusters;
	in_str >> num_clusters;
	num_clusters /= N;

	std::cout << "Number of clusters in file: " << num_clusters << "\n";

	//create an array of HC objects
	HCC* hc = new HCC(num_clusters, N, maxT);

	double x, y;
	double temp;
	int particle = 0;
	int cluster = 0;
	int timestep = 0;


	//fill the hcc with clusters
	while (in_str >> x) {
		//x is stored already, store y
		in_str >> y;

		//determine where in the cluster array these coordinates goes
		int placement = DIMENSION * particle + timestep * N * DIMENSION;
		(*hc)[cluster].clusters[placement] = x;
		(*hc)[cluster].clusters[placement+1] = y;

		//update indices
		particle++;
		if (particle == N) {              //increment cluster by 1 every N particles
			cluster++; particle = 0;
			if (cluster == num_clusters) {  //increment timestep by 1 after all clusters
				timestep++; cluster = 0;
				in_str >> temp;               //extra value when timestep changes
				if (timestep == maxT) {       //break if run out of storage
					break;                      
				}
			}
		}

		//cycle through the useless data
		for (int i = 0; i < 5; i++) {
			in_str >> temp;
		}

		//std::cout << "Time step " << timestep << "\n";

	}

	//close the input file
	in_str.close();

	//return the collection
	return hc;
}

// ********************************************************************** //
// ********************* Transition Detection *************************** //
// ********************************************************************** //


void getCoordinates(HydroCluster& hc, double* X, int N, int time) {
	//place the coordinates of the cluster at given time into X

	int first_element = time * N * DIMENSION; 
	for (int i = 0; i < N * DIMENSION; i++) {
		X[i] = hc.clusters[first_element + i];
	}
}

void getAdj(double* X, int N, int* M, double cutoff) {
	//get the adjacnecy matrix from a cluster X - specify cutoff

	//init a particle array to 0, then fill it from X
	double* particles = new double[DIMENSION*N]; double* Z = new double[DIMENSION];
	for (int p = 0; p < DIMENSION*N; p++) {
		particles[p] = 0;
	}
	c2p(X, particles, N);

	/*
	for (int q = 0; q < 2*N; q++) {
		std::cout << particles[q] << "\n";
	}
	*/

	//compute the distance between particles, construct adj matrix
	int i, j; 
	for (int index = 0; index < N*N; index++) {
		index2ij(index, N, i, j);
		if (j > i) {
			double d = euDist(particles, i, j, N, Z);
			if (d < cutoff) {
				M[index] = 1; M[toIndex(j,i,N)] = 1;
			}
		}
	}

	delete []particles; delete []Z;
}

void checkState(int N, double* X, int state, Database* db,
							  bool& reset, int& new_state) {
	//check if a new state has been reached

	//declare reset false, only true if the sample is invalid
	reset = false;

	//get adjacency matrix of the current state
	int* M = new int[N*N]; for (int i = 0; i < N*N; i++) M[i]=0;
	getAdj(X, N, M, HYDRO_CUT);

	//compare old and new state, check if same
	int* old = new int[N*N]; for (int i = 0; i < N*N; i++) old[i]=0;
	int old_bonds = (*db)[state].getBonds();
	extractAM(N, state, old, db);
	int success = checkSame(old, M, N);

	int B = 0;
	for (int i = 0; i < N*N; i++) {
		B += M[i];
	}
	B /= 2;
	std::cout << B << "\n";

	if (success) {//same matrix
		delete []M; delete []old;
		return;
	}
	else {//not the same, find matrix in database
		success = mcm::findMatrix(M, old, old_bonds, N, db,  reset,  new_state);
		if (!success) {
			//may output an unphysical state. refine with newton, check again
			refine(N, X, M);
			success = mcm::findMatrix(M, old, old_bonds, N, db,  reset,  new_state);
			if (!success) { //state still not found after refine, ignore this sample
				reset = true; 
				printf("State not found in database, B = %d \n\n\n\n\n", B);
			}
		}
		delete []M; delete []old;
	}
}

//maybe parralelize this later?
void determineTransitions(HCC* hc, Database* db) {
	//takes a collection of cluster trajectories, determines the folding pathway and
	//times, constructs an adjacency matrix and time distribution

	//get info from the db and HCC
	int N = db->getN();             //number of particles
	int ns = db->getNumStates();    //number of states
	int nc = hc->getNumClusters();  //number of clusters
	int maxT = hc->getMaxT();       //max timestep for clusters evo

	//declare all the things we need
	double* X = new double[N*DIMENSION];    //store a configuration
	double* Xold = new double[N*DIMENSION]; //store previous config
	int timer = 0;
	std::vector<int>* transition_times = new std::vector<int>[ns*ns];
	int old_state;
	int new_state;
	bool reset;
	int bonds;

	//loop over each cluster in HHC
	for (int i = 8; i < 9; i++) {
		//create a reference to the i-th collection of clusters
		HydroCluster& currentCluster = (*hc)[i];
		printf("Now analyzing cluster %d\n", i);

		//put the 0-th cluster in temp - get details
		getCoordinates(currentCluster, Xold, N, 0);
		old_state = currentCluster.cNum; 
		new_state = currentCluster.cNum; 
		//for (int k = 0; k < 12; k++) std::cout << Xold[k] << "\n";

		//loop over all the timesteps - detect transitions
		for (int time = 1; time < maxT; time++) {

			//get the next set of coordinates - incremenent timer
			getCoordinates(currentCluster, X, N, time);
			if (X[0] == 0) {
				std::cout << time << "\n";
				break;
			}
			timer++;

			//check if state changed
			checkState(N, X, old_state, db, reset, new_state);


			if (reset) {
				std::cout << "something bad happened\n";
				timer--;
			}

			if (new_state != old_state) { //new state reached, update data
				printf("Transition %d to %d in time %d\n", old_state, new_state, timer);
				transition_times[toIndex(old_state, new_state, ns)].push_back(timer);
				old_state = new_state;
				timer = 0;

				//check if ground state has been reached
				bonds = (*db)[new_state].getBonds();
				if (bonds >= 2*N-3) {
					break;
				}
			}
			
			//copy new state to old state
			if (!reset) {
				memcpy(Xold, X, DIMENSION*N*sizeof(double)); // copy X to Xold
			}

		}
	}

	//free the memory
	delete []X; delete []Xold; delete []transition_times;


}




// ********************************************************************** //
// ********************* Testing Area  ********************************** //
// ********************************************************************** //


void testExtract(HCC* hc) {
	//test if the extraction is working

	int test_cluster = 0;
	int time = 6;

	for (int i = 0; i < 12; i++) {
		std::cout << (*hc)[test_cluster].clusters[time*12+i] << "\n";
	}
}










}