#include "database.h"
#include "nauty.h"
#include "hydro.h"
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include "../defines.h"
namespace bd {

HydroCluster::HydroCluster() {
	//constructor

	cNum = 0;
	t_last = 0;
	double* clusters = new double[N*DIMENSION*maxT];
}

HydroCluster::~HydroCluster() {
	//deconstructor

	delete []clusters;
}

HCC::HCC(int num_clusters_, int N_, int maxT_) {
	//constructor

	num_clusters = num_clusters_;
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

	std::cout << "num clusters " << num_clusters << "\n";

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

//maybe parralelize this later?
void determineTransitions(HCC* hc, Database* db) {
	//takes a collection of cluster trajectories, determines the folding pathway and
	//times, constructs and adjacency matrix and time distribution



}




void testExtract(HCC* hc) {
	//test if the extraction is working

	int test_cluster = 0;

	for (int i = 1100; i < 1120; i++) {
		std::cout << (*hc)[test_cluster].clusters[12*i] << (*hc)[test_cluster].clusters[12*i+1] <<
		(*hc)[test_cluster].clusters[12*i+2] << (*hc)[test_cluster].clusters[12*i+3] << 
		(*hc)[test_cluster].clusters[12*i+4] << (*hc)[test_cluster].clusters[12*i+5] << "\n\n";
	}
}










}