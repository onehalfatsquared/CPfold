#include "design.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <map>
#include <eigen3/Eigen/Dense>
#include "point.h"
#include "pair.h"
#include "database.h"
#include "bDynamics.h"
#include "tpt.h"
#include "nauty.h"
#include "graph.h"
#include "graphviz.h"
#include "../defines.h"
namespace bd{




int setTypes(int N, int* particleTypes, int IC) {
	//set which particles are which type

	int numTypes;

	if (IC == 0) { //all the same type
		for (int i = 0; i < N; i++) {
			particleTypes[i] = 0;
		}
		numTypes = 1;
	}
	else if (IC == 1) { //alternating type a and b
		for (int i = 0; i < N; i++) {
			particleTypes[i/2] = 0;
		}
		for (int i = 0; i < N; i++) {
			particleTypes[i/2+1] = 1;
		}
		numTypes = 2;
	}

	return numTypes;
}

int readDesignFile(int N, int* particleTypes) {
	//read the file to get particle identities

	int max = -1; //use to determine number of types
	int count = 0; //use to determine if input mathces the system

	//file location
	std::string filename = "input/design/particles.txt";

	//open the file
	std::ifstream in_str(filename);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		return -1;
	}

	//store the entries
	int id; 
	while (in_str >> id) {
		//increment the counter for the particle number (1 based indexing)
		count++; 

		//if the particle number doesnt exceed the total, store it
		if (count <= N) {
			particleTypes[count-1] = id;
		}

		//get the max of the particle ids
		if (id > max) {
			max = id;
		}
	}

	//check if the number of entries in file matches the db we recieved
	if (count != N) {
		fprintf(stderr, "The supplied design file is incompatible with the database file.");
		return -1;
	}

	//return number of particle types
	return max+1;
}

double getStickyProduct(int N, int state,  Database* db, int* particleTypes, 
												std::map<std::pair<int,int>,double> kappa) {
	//return the product of sticky param ^ num bonds for each type

	double stickyProd = 1;

	//loop over the adjacency matrix, determine bond types, get factors
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			if ((*db)[state].isInteracting(i,j,N)) {
				int p1 = particleTypes[i]; int p2 = particleTypes[j];
				stickyProd *= kappa[{p1,p2}];
			}
		}
	}

	return stickyProd;
}

void reweight(int N, int num_states, Database* db, int* particleTypes, double* eq,
						  std::map<std::pair<int,int>,double> kappa) {
	//perform the re-weighting of the eq measure for new sticky params

	double kap0 = KAP;                 //sticky parameter for initial measurement
	double beta = BETA;                //inverse temp

	double Z = 0;                     //normalizing constant for eq

	//loop over states, get the equilibrium measure entry
	for (int i = 0; i < num_states; i++) {
		//get initial eq prob and number of bonds
		int b = (*db)[i].getBonds();
		double prob = (*db)[i].getFrequency();

		//get the denominator of reweight factor - KAP^b_i
		double denom = pow(KAP, b*BETA);

		//get the numerator of the re-weight factor
		double num = pow(getStickyProduct(N, i, db, particleTypes, kappa),BETA); 

		//compute the new probability
		double new_prob = prob * num / denom;

		//debug line for corect reweight
		//printf("Test: Old %f, New %f\n", prob, new_prob);

		//increment Z and add to array
		Z += new_prob;
		eq[i] = new_prob;
	}

	//re-normalize
	for (int i = 0; i < num_states; i++) eq[i] /= Z;
}

void initKappaVals(int numInteractions, double* kappaVals) {
	//initially set the valyes to all 1 for each interaction

	for (int i = 0; i < numInteractions; i++) {
		kappaVals[i] = 1.0;
	}
}

void makeKappaMap(int numTypes, double* kappaVals, 
									std::map<std::pair<int,int>,double>& kappa) {
	//fill the map with values in kappaVals

	int counter = 0; 

	for (int i = 0; i < numTypes; i++) {
		for (int j = i; j < numTypes; j++) {
			kappa[{i,j}] = kappaVals[counter];
			kappa[{j,i}] = kappaVals[counter];
			counter++;
		}
	}

}

void getGroundStates(int N, Database* db, std::vector<int>& gs) {
	//return a vector with all the ground states in it

	int ns = db->getNumStates(); 
	for (int i = 0; i < ns; i++) {
		int b = (*db)[i].getBonds();
		if (b == 2*N-3 || b == 2*N-2) {
			gs.push_back(i);
		}
	}
}



void constructSurface(int N, Database* db, int initial, int target, bool useFile) {
	//plot the surface of hitting probabilities for initial to target states

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter map
	double* kappaVals = new double[numInteractions];
	initKappaVals(numInteractions, kappaVals);
	std::map<std::pair<int,int>,double> kappa;
	makeKappaMap(numTypes, kappaVals, kappa);

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);

	//get all the ground states
	std::vector<int> ground;
	getGroundStates(N, db, ground);

	//create the reweighted equilibrium measure
	double* eq = new double[num_states];
	reweight(N, num_states, db, particleTypes, eq, kappa); 

	//do hitting probability calculation




	//free memory
	delete []particleTypes; delete []kappaVals; delete []eq;
}




}