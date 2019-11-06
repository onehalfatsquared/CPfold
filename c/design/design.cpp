#include "design.h"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <map>
#include <deque>
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

void allPerms(int N, std::deque<std::string>& perms) {
	//generate all permutations of N 0s and 1s

	std::string p0 = "0"; std::string p1 = "1";
	perms.push_back(p0); perms.push_back(p1);
	int elements = 2;
	int final_size = pow(2,N); 
	while (elements < final_size) {
		std::string f = perms.front(); perms.pop_front();
		perms.push_back(f+p0); perms.push_back(f+p1);
		elements += 1;
	}
}

void setTypes(int N, int* particleTypes, std::deque<std::string> perms, int perm) {
	//select a permutation from the perms vector and set particle types

	std::string p = perms[perm];
	for (int i = 0; i < N; i++) {
		particleTypes[i] = p[i]-'0';
		std::cout << particleTypes[i] << "\n";
	}
}


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
		kappaVals[i] = 3.0;
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

double getHitProbability(int num_states, int initial, std::vector<int> targets, double* U) {
	//get the hitting probability for the group of states in targets

	double prob = 0; 
	for (int i = 0; i < targets.size(); i++) {
		prob += U[toIndex(initial, targets[i], num_states)];
	}

	return prob;
}

double solveAbsorbProblem(int initial, double* kappaVals, Database* db, int* particleTypes, 
								 					double* Tconst, std::vector<int> ground, std::vector<int> targets) {
	//set up and solve the absorbtion probability problem

	//get parameters from database
	int N = db->getN();
	int num_states = db->getNumStates();

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* T = new double[num_states*num_states]; //rate matrix
	double* P = new double[num_states*num_states]; //probability transition matrix
	double* U = new double[num_states*num_states]; //hitting probability matrix
	double* eq = new double[num_states];           //equilibrium measure

	//declare the kappa map
	std::map<std::pair<int,int>,double> kappa;

	//compute the hitting probability at the current kappa vals
	//make the map
	makeKappaMap(2, kappaVals, kappa);
	//do rewieght
	reweight(N, num_states, db, particleTypes, eq, kappa);
	//reset transition matrices
	for (int i = 0; i < num_states*num_states; i++) {
		P[i] = U[i] = 0;
	}
	//copy Tconst into T
	std::copy(Tconst, Tconst+num_states*num_states, T);
	//fill in transposed entries such that T satisfies detailed balance
	satisfyDB(T, num_states, db, eq);
	//fill in diagonal with negative sum of row entries
	fillDiag(T, num_states);
	//use the filled rate matrix to compute probability matrix
	createProbabilityMatrix(T, num_states, P);
	//solve for hitting probabilities to endStates states
	computeHittingProbability(P, num_states, ground, U);
	//get hit prob to targets
	double p = getHitProbability(num_states, initial, targets, U);

	delete []T; delete []P; delete []U; delete []eq;

	return p;
}

double computeGrad(int initial, int numInteractions, double* kappaVals, Database* db, 
								 int* particleTypes, double* Tconst, std::vector<int> ground, 
								 std::vector<int> targets, double* g) {
	//compute the gradient of the hitting probability in kappa
	//returns the hitting probability at current point

	//declare parameters
	double h = 1e-5; //step size for finite differences

	//get the hitting probability at the current kappa vals
	double f0 = solveAbsorbProblem(initial, kappaVals, db, particleTypes, Tconst, ground, targets);

	//std::cout << "Hitting Prob = " << f0 << "\n";

	//take a step h in each direction and compute new hitting prob
	for (int k = 0; k < numInteractions; k++) {
		if (kappaVals[k] < 0.011) {
			g[k] = 0;
		}
		else {
			kappaVals[k] += h;
			double f1 = solveAbsorbProblem(initial, kappaVals, db, particleTypes, Tconst, ground, targets);
			kappaVals[k] -= h;

			g[k] = (f1 - f0) / h;
		}
	}

	return f0;
}

void hittingProbMaxTOYperms(int N, Database* db, int initial, int target) {
	//use optimization scheme to get max hitting probability for target in toy model
	//does optimization over all particle perms and outputs max probability

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 2000; double tol = 2e-3; double step = 0.9;

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes = 2;
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//init a gradient array
	double* g = new double[numInteractions];

	//loop over all permutations
	std::deque<std::string> perms; allPerms(N,perms);
	int num_perms = perms.size();
	double* permProb = new double[num_perms];

	for (int p = 0; p < num_perms; p++) {
		//get the permutation
		setTypes(N, particleTypes, perms, p);
		initKappaVals(numInteractions, kappaVals);
		double hit;

		//optimization - steepest ascent
		for (int it = 0; it < max_its; it++) {
			//compute the gradient of kappa
			hit = computeGrad(initial, numInteractions, kappaVals, db, particleTypes, 
												Tconst,  ground, targets, g);

			//update kappa - set to 0.1 if it goes negative
			for (int k = 0; k < numInteractions; k++) {
				kappaVals[k] += step*g[k];
				if (kappaVals[k] < 0) {
					kappaVals[k] = 0.01;
				}
			}

			//check for gradient tolerance
			double res = 0;
			for (int k = 0; k < numInteractions; k++) res += abs(g[k]);
			std::cout << "Iteration " << it << ", Res: " << res << "\n";
			if (res < tol) {
				break;
			}
		}

		permProb[p] = hit;
	}

	//print out the hitting probability with permutation
	for (int p = 0; p < num_perms; p++) {
		printf("Hitting Prob: %f, ID: %s\n", permProb[p], perms[p].c_str());
	}

	//std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []Tconst; delete []g; delete []permProb;
}

void hittingProbMaxTOY(int N, Database* db, int initial, int target, bool useFile) {
	//use optimization scheme to get max hitting probability for target in toy model
	//uses the file for the permutation

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 2000; double tol = 2e-3; double step = 0.9;

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

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//init a gradient array
	double* g = new double[numInteractions];

	//get the permutation
	initKappaVals(numInteractions, kappaVals);
	double hit;

	//optimization - steepest ascent
	for (int it = 0; it < max_its; it++) {
		//compute the gradient of kappa
		hit = computeGrad(initial, numInteractions, kappaVals, db, particleTypes, 
											Tconst,  ground, targets, g);

		//update kappa - set to 0.1 if it goes negative
		for (int k = 0; k < numInteractions; k++) {
			kappaVals[k] += step*g[k];
			if (kappaVals[k] < 0) {
				kappaVals[k] = 0.01;
			}
		}

		//check for gradient tolerance
		double res = 0;
		for (int k = 0; k < numInteractions; k++) res += abs(g[k]);
		std::cout << "Iteration " << it << ", Res: " << res << "\n";
		if (res < tol) {
			break;
		}
	}

	std::cout << "Maximum Probability: " << hit << "\n";

	std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []Tconst; delete []g; 
}



void constructSurfaceTOY(int N, Database* db, int initial, int target, bool useFile) {
	//plot the surface of hitting probabilities for initial to target states
	//toy model where 3 sticky parameters is hard coded - 2 unique

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

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];
	initKappaVals(numInteractions, kappaVals);

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* T = new double[num_states*num_states]; //rate matrix
	double* P = new double[num_states*num_states]; //probability transition matrix
	double* U = new double[num_states*num_states]; //hitting probability matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries
	double* eq = new double[num_states];           //equilibrium measure

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//declare the kappa mapping
	std::map<std::pair<int,int>,double> kappa;

	//declare outfile
	std::ofstream ofile;
	ofile.open("hittingProbSurface.txt");

	//do hitting probability calculation
	int M = 60; //num points in each dimension
	for (int x = 0; x < M; x++) {
		for (int y = 0; y < M; y++) {
			//set kappa
			kappaVals[0] = 0.5+ (double)x / 4; 
			kappaVals[1] = 0.5+ (double)y / 4; 
			kappaVals[2] = 0.5+ (double)x / 4; 

			//make the map
			makeKappaMap(2, kappaVals, kappa);
			//std::cout << kappa[{0,0}] << ' ' << kappa[{1,0}] << ' ' << kappa[{1,1}] << "\n";
			//do rewieght
			reweight(N, num_states, db, particleTypes, eq, kappa);
			//reset transition matrices
			for (int i = 0; i < num_states*num_states; i++) {
				P[i] = U[i] = 0;
			}
			//copy Tconst into T
			std::copy(Tconst, Tconst+num_states*num_states, T);
			//fill in transposed entries such that T satisfies detailed balance
			satisfyDB(T, num_states, db, eq);
			//fill in diagonal with negative sum of row entries
			fillDiag(T, num_states);
			//use the filled rate matrix to compute probability matrix
			createProbabilityMatrix(T, num_states, P);
			//solve for hitting probabilities to endStates states
			computeHittingProbability(P, num_states, ground, U);
			//get hit prob to targets
			double p = getHitProbability(num_states, initial, targets, U);

			//output to file
			ofile << kappaVals[0] << ' ' << kappaVals[1] << ' ' << p << "\n";

		}
		std::cout << x << "\n";
	}




	//close file
	ofile.close();

	//free memory
	delete []particleTypes; delete []kappaVals; delete []eq;
	delete []T; delete []P;  delete []U; delete []Tconst;
}




}