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
#include "adjacency.h"
#include "database.h"
#include "bDynamics.h"
#include "tpt.h"
#include "nauty.h"
#include "graph.h"
#include "graphviz.h"
#include "../defines.h"
namespace bd{

/****************************************************/
/*********** General Functions    *******************/
/****************************************************/

void getBondTypes(int N, int* particleTypes, Database* db, std::vector<int> targets) {
	//for each state in targets, check how many of each bond type each state has
	//hard-coded for 2 types

	for (int index = 0; index < targets.size(); index++) {
		int AA = 0; int AB = 0; int BB = 0;
		int state = targets[index];
		for (int i = 0; i < N; i++) {
			for (int j = i+2; j < N; j++) {
				if ( (*db)[state].isInteracting(i,j,N)) {
					int p1 = particleTypes[i];
					int p2 = particleTypes[j];
					if (p1 == 0 && p2 == 0) {
						AA++;
					}
					else if (p1 == 1 && p2 == 1) {
						BB++;
					}
					else {
						AB++;
					}
				}
			}
		}

		printf("State: %d, AA: %d, AB %d, BB %d\n", state, AA, AB, BB);
	}
}

void checkPerm(std::string word, std::deque<std::string>& perms) {
	//check if reflected and swapped word is in perms, if not add it

	//check reflections
	std::string rword = word; std::reverse(rword.begin(), rword.end());
	std::deque<std::string>::iterator itR = std::find(perms.begin(), perms.end(), rword);
	if (itR != perms.end()) {
		return;
	}

	//check swaps
	std::string swapWord = word; 
	for (int i = 0; i < word.length(); i++) {
		if (word[i] == '0') swapWord[i] = '1';
		if (word[i] == '1') swapWord[i] = '0';
	}
	std::deque<std::string>::iterator itS = std::find(perms.begin(), perms.end(), swapWord);
	if (itS != perms.end()) {
		return;
	}

	//check reflections of swaps
	std::string rword2 = swapWord; std::reverse(rword2.begin(), rword2.end());
	std::deque<std::string>::iterator itR2 = std::find(perms.begin(), perms.end(), rword2);
	if (itR2 != perms.end()) {
		return;
	}

	perms.push_back(word);
}

void allPerms(int N, std::deque<std::string>& perms) {
	//generate all perms of N 0s and 1s

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

void distinctPerms(int N, std::deque<std::string>& perms) {
	//generate all permutations of N 0s and 1s w/o reflections

	//first generate all (N-1) length possibilities, fill preliminary deque
	std::string p0 = "0"; std::string p1 = "1";
	std::deque<std::string> prelim;
	prelim.push_back(p0); prelim.push_back(p1);
	int elements = 2;
	int final_size = pow(2,N-1); 
	while (elements < final_size) {
		std::string f = prelim.front(); prelim.pop_front();
		prelim.push_back(f+p0); prelim.push_back(f+p1);
		elements += 1;
	}

	//next, add 0 and 1 to each (N-1) length permuation, only add it to perms
	//if its reflection is not already in it
	for (int i = 0; i < prelim.size(); i++) {
		std::string f = prelim[i];
		checkPerm(f+p0, perms); checkPerm(f+p1, perms);
	}

	printf("Found %ld distinct permutations\n", perms.size());

}

void setTypes(int N, int* particleTypes, std::deque<std::string> perms, int perm) {
	//select a permutation from the perms vector and set particle types

	std::string p = perms[perm];
	for (int i = 0; i < N; i++) {
		particleTypes[i] = p[i]-'0';
		//std::cout << particleTypes[i] << "\n";
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
		kappaVals[i] = 3;
	}
}

void checkPositive(int numInteractions, double* kappaVals) {
	//check if each kappa is positive, if not set to 0.01

	for (int i = 0; i < numInteractions; i++) {
		if (kappaVals[i] < 0) {
			kappaVals[i] = 0.01;
		}
	}
}

void applyMax(int numInteractions, double* kappaVals, double M) {
	//check if each kappa is below max, M. set to M-0.01 if not

	for (int i = 0; i < numInteractions; i++) {
		if (kappaVals[i] > M) {
			kappaVals[i] = M-0.01;
		}
	}
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

/****************************************************/
/******** Hitting Probability Optimization   ********/
/****************************************************/

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



double computeGradHP(int initial, int numInteractions, double* kappaVals, Database* db, 
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
	std::deque<std::string> perms; 
	//allPerms(N,perms);
	distinctPerms(N, perms);
	int num_perms = perms.size();
	double* permProb = new double[num_perms];

	for (int p = 0; p < num_perms; p++) {
		//get the permutation
		setTypes(N, particleTypes, perms, p);
		initKappaVals(numInteractions, kappaVals);
		double hit;

		printf("Testing permutation %d of %d\n", p+1, num_perms);

		//optimization - steepest ascent
		for (int it = 0; it < max_its; it++) {
			//compute the gradient of kappa
			hit = computeGradHP(initial, numInteractions, kappaVals, db, particleTypes, 
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
		hit = computeGradHP(initial, numInteractions, kappaVals, db, particleTypes, 
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

	double R = getRate(initial, kappaVals, db, particleTypes, Tconst, targets);
	printf("The average transition rate with this kappa is %f\n", R);


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []Tconst; delete []g; 
}



/****************************************************/
/******** Equilinrium Probability Optimization   ****/
/****************************************************/


double getEqProb(int initial, double* kappaVals, Database* db, 
								 int* particleTypes, std::vector<int> targets) {
	//do re-weighting and get the eqiulibrium probability sum of all target states

	//get parameters from database
	int N = db->getN();
	int num_states = db->getNumStates();

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* eq = new double[num_states];           //equilibrium measure

	//declare the kappa map
	std::map<std::pair<int,int>,double> kappa;

	//compute the eq probability at the current kappa vals
	//make the map
	makeKappaMap(2, kappaVals, kappa);
	//do rewieght
	reweight(N, num_states, db, particleTypes, eq, kappa);

	//get the new equilibrium probability of target state
	double prob = 0;
	for (int i = 0;  i< targets.size(); i++) {
		prob += eq[targets[i]];
	}

	//free memory
	delete []eq;

	return prob;
}

double computeGradEQ(int initial, int numInteractions, double* kappaVals, Database* db, 
								 int* particleTypes, std::vector<int> targets, double* g) {
	//compute the gradient of the eqiluibrium probability in kappa
	//returns the eq probability at current point

	//declare parameters
	double h = 1e-5; //step size for finite differences

	//get the hitting probability at the current kappa vals
	double f0 = getEqProb(initial, kappaVals, db, particleTypes, targets);

	//std::cout << "Hitting Prob = " << f0 << "\n";

	//take a step h in each direction and compute new hitting prob
	for (int k = 0; k < numInteractions; k++) {
		if (kappaVals[k] < 0.011) {
			g[k] = 0;
		}
		else if (kappaVals[k] > 100) {
			g[k] = 0;
		}
		else {
			kappaVals[k] += h;
			double f1 = getEqProb(initial, kappaVals, db, particleTypes, targets);
			kappaVals[k] -= h;

			g[k] = (f1 - f0) / h;
		}
	}

	return f0;
}

void eqProbMaxTOY(int N, Database* db, int initial, int target, bool useFile) {
	//use optimization scheme to get max equilibrium probability for target in toy model
	//uses the file for the permutation

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 120000; double tol = 1e-6; double step = 0.9;

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

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	getBondTypes(N, particleTypes, db, targets);

	return;

	//init a gradient array
	double* g = new double[numInteractions];

	//get the permutation
	initKappaVals(numInteractions, kappaVals);
	//readKappaFile(numInteractions, kappaVals);
	double eq;

	//optimization - steepest ascent
	for (int it = 0; it < max_its; it++) {
		//compute the gradient of kappa
		eq = computeGradEQ(initial, numInteractions, kappaVals, db, particleTypes, 
											targets, g);

		//update kappa - set to 0.1 if it goes negative
		for (int k = 0; k < numInteractions; k++) {
			kappaVals[k] += step*g[k];
			if (kappaVals[k] < 0) {
				kappaVals[k] = 0.01;
			}
			else if (kappaVals[k] > 100) {
				kappaVals[k] = 99.9;
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

	std::cout << "Maximum Eq Probability: " << eq << "\n";

	std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []g; 
}

void eqProbMaxTOYperms(int N, Database* db, int initial, int target) {
	//use optimization scheme to get max hitting probability for target in toy model
	//does optimization over all particle perms and outputs max probability

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 60000; double tol = 1e-6; double step = 0.9;

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes = 2;
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//init a gradient array
	double* g = new double[numInteractions];

	//loop over all permutations
	std::deque<std::string> perms; 
	//allPerms(N, perms);
	distinctPerms(N,perms);
	int num_perms = perms.size();
	double* permProb = new double[num_perms];

	for (int p = 0; p < num_perms; p++) {
		//get the permutation
		setTypes(N, particleTypes, perms, p);
		initKappaVals(numInteractions, kappaVals);
		double eq;

		printf("Testing permutation %d of %d\n", p+1, num_perms);

		//optimization - steepest ascent
		for (int it = 0; it < max_its; it++) {
			//compute the gradient of kappa
			eq = computeGradEQ(initial, numInteractions, kappaVals, db, particleTypes, 
												targets, g);

			//update kappa - set to 0.1 if it goes negative
			for (int k = 0; k < numInteractions; k++) {
				kappaVals[k] += step*g[k];
				if (kappaVals[k] < 0) {
					kappaVals[k] = 0.01;
				}
				else if (kappaVals[k] > 100) {
				kappaVals[k] = 99.9;
				}
			}

			//check for gradient tolerance
			double res = 0;
			for (int k = 0; k < numInteractions; k++) res += abs(g[k]);
			//std::cout << "Iteration " << it << ", Res: " << res << "\n";
			if (res < tol) {
				break;
			}
		}

		permProb[p] = eq;
	}

	//print out the hitting probability with permutation
	for (int p = 0; p < num_perms; p++) {
		printf("Eq Prob: %f, ID: %s\n", permProb[p], perms[p].c_str());
	}

	//std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []g; delete []permProb;
}


/****************************************************/
/******** Transition Rate  Optimization   ***********/
/****************************************************/

double getRate(int initial, double* kappaVals, Database* db, int* particleTypes, 
							 double* Tconst, std::vector<int> targets) {
	//do re-weighting and get the transition rate to the target state

	//get parameters from database
	int N = db->getN();
	int num_states = db->getNumStates();

	//initialize the transition rate matrix
	double* T = new double[num_states*num_states];
	//copy Tconst into T
	std::copy(Tconst, Tconst+num_states*num_states, T);

	//init array for equilibrium distribution and compute it
	double* eq = new double[num_states];
	//declare and fill kappa map
	std::map<std::pair<int,int>,double> kappa;
	makeKappaMap(2, kappaVals, kappa);   
	//do reweight to fill eq
	reweight(N, num_states, db, particleTypes, eq, kappa);

	//fill in transposed entries such that T satisfies detailed balance
	satisfyDB(T, num_states, db, eq);

	//fill in diagonal with negative sum of row entries
	fillDiag(T, num_states);

	//solve for the committor
	//initialize committor in q
	double* q = new double[num_states];
	for (int i = 0; i < num_states; i++) q[i] = 0;
	
	//solve dirichlet problem for committor, q 
	computeCommittor(q, T, num_states, initial, targets);

	//get the average transition rate
	double R = computeTransitionRate(num_states, q, T, eq);

	//free memory 
	delete []T; delete []eq; delete []q;

	return R;
}

double computeGradRate(int initial, int numInteractions, double* kappaVals, Database* db, 
								 int* particleTypes, double* Tconst, std::vector<int> targets, double* g) {
	//compute the gradient of the avg transition rate in kappa
	//returns the avg transition rate at current point

	//declare parameters
	double h = 1e-5; //step size for finite differences

	//get the hitting probability at the current kappa vals
	double f0 = getRate(initial, kappaVals, db, particleTypes, Tconst, targets);

	//std::cout << "Hitting Prob = " << f0 << "\n";

	//take a step h in each direction and compute new hitting prob
	for (int k = 0; k < numInteractions; k++) {
		if (kappaVals[k] < 0.011) {
			g[k] = 0;
		}
		else {
			kappaVals[k] += h;
			double f1 = getRate(initial, kappaVals, db, particleTypes, Tconst, targets);
			kappaVals[k] -= h;

			g[k] = (f1 - f0) / h;
		}
	}

	return f0;
}



void rateMaxTOY(int N, Database* db, int initial, int target, bool useFile) {
	//use optimization scheme to get max rate for target in toy model
	//uses the file for the permutation

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 5000; double tol = 1e-10; double step = 0.9;

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
	//initKappaVals(numInteractions, kappaVals);
	readKappaFile(numInteractions, kappaVals);
	double R;

	//optimization - steepest ascent
	for (int it = 0; it < max_its; it++) {
		//compute the gradient of kappa
		R = computeGradRate(initial, numInteractions, kappaVals, db, particleTypes, 
											Tconst, targets, g);

		//update kappa - set to 0.1 if it goes negative
		for (int k = 0; k < numInteractions; k++) {
			kappaVals[k] += step*g[k];
		}
		checkPositive(numInteractions, kappaVals);
		applyMax(numInteractions, kappaVals, 100);

		//check for gradient tolerance
		double res = 0;
		for (int k = 0; k < numInteractions; k++) res += abs(g[k]);
		std::cout << "Iteration " << it << ", Res: " << res << "\n";
		if (res < tol) {
			break;
		}
	}

	std::cout << "Maximum Rate: " << R << "\n";

	std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";

	double hp = solveAbsorbProblem(initial, kappaVals, db, particleTypes, Tconst, 
																 ground, targets);
	printf("The hitting probability at this kappa is %f\n", hp);


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []g; delete []Tconst;
}

void rateMaxTOYperms(int N, Database* db, int initial, int target) {
	//use optimization scheme to get max rate for target in toy model
	//uses the file for the permutation

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 5000; double tol = 1e-10; double step = 0.9;

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
	std::deque<std::string> perms; 
	//allPerms(N,perms);
	distinctPerms(N, perms);
	int num_perms = perms.size();
	double* permRate = new double[num_perms];

	for (int p = 0; p < num_perms; p++) {
		//get the permutation
		setTypes(N, particleTypes, perms, p);
		initKappaVals(numInteractions, kappaVals);
		double R;

		printf("Testing permutation %d of %d\n", p+1, num_perms);

		//optimization - steepest ascent
		for (int it = 0; it < max_its; it++) {
			//compute the gradient of kappa
			R = computeGradRate(initial, numInteractions, kappaVals, db, particleTypes, 
												Tconst, targets, g);

			//update kappa - set to 0.1 if it goes negative
			for (int k = 0; k < numInteractions; k++) {
				kappaVals[k] += step*g[k];
			}
			checkPositive(numInteractions, kappaVals);
			applyMax(numInteractions, kappaVals, 100);

			//check for gradient tolerance
			double res = 0;
			for (int k = 0; k < numInteractions; k++) res += abs(g[k]);
			//std::cout << "Iteration " << it << ", Res: " << res << "\n";
			if (res < tol) {
				break;
			}
		}

		permRate[p] = R;
	}

	//print out the hitting probability with permutation
	for (int p = 0; p < num_perms; p++) {
		printf("Average Rate: %f, ID: %s\n", permRate[p], perms[p].c_str());
	}

	//std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []g; delete []permRate; delete []Tconst;
}


/****************************************************/
/*** Hit Prob Optimization w/ rate constraint   *****/
/****************************************************/

void applyRateConstraint(double R, double c, double r, int numInteractions, 
												 double* gH, double* gR) {
	//determine if the constraint needs to be applied and return appropriate gradient

	if (R < c) {//constraint active
		for (int i = 0; i < numInteractions; i++) gH[i] += r*gR[i];
	}

}


double lineSearch(int initial, int numInteractions, double* kappaVals, Database* db, int* particleTypes, 
									double* Tconst, std::vector<int> ground, std::vector<int> targets, 
									double c, double r, double H, double R, double* g, double& step) {
	//perform a line search to get a step size for the objective fn
	//step is passed by reference, returns difference in obj fn

	//define initial step size and reducing factor
	double step0 = 2; double alpha = 0.5;

	//see if the constraint needs to be applied
	int apply = 0;
	if (R < c) {
		apply = 1;
	}

	//evaluate objective fn
	double obj = H + r*apply*(R-c);

	//printf("Original Objective: %f\n", obj);

	//declare array for test kappaVals
	double* testKappa = new double[numInteractions];


	//perform line searches
	double objNew; double delta = -1; 
	int maxReductions = 20; 
	for (int i = 0; i < maxReductions; i++) {
		//get the new kappa values
		step0 *= alpha;
		for (int k = 0; k < numInteractions; k++) testKappa[k] = kappaVals[k] + g[k]*step0;
		checkPositive(numInteractions, testKappa);

		//evaluate the objective
		double Hnew = solveAbsorbProblem(initial, testKappa, db, particleTypes, Tconst, 
																	 ground, targets);
		double Rnew = getRate(initial, testKappa, db, particleTypes, Tconst, targets);
		int applyNew = 0;
		if (Rnew < c) {
			applyNew = 1;
		}
		objNew = Hnew + r*applyNew*(Rnew-c);
		//printf("New obj: %f\n", objNew);

		delta = objNew - obj;
		if (delta > 0) {
			break;
		}
	}

	step = step0;


	//free memory
	delete []testKappa;

	return delta;

}

void hittingProbMaxTOYc(int N, Database* db, int initial, int target, bool useFile) {
	//use optimization scheme to get max hitting probability for target in toy model
	//has constraint for transition rate
	//uses the file for the permutation

	//get database info
	int num_states = db->getNumStates(); 

	//set iteration settings
	int max_its = 2000; double objTol = 1e-5; 
	double c = 0.12; //rate must be > c (target dependent)
	double r = 200;    //cost 

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
	double* gH = new double[numInteractions];
	double* gR = new double[numInteractions];

	//get the permutation
	//initKappaVals(numInteractions, kappaVals);
	readKappaFile(numInteractions, kappaVals);
	double hit;
	double R; 

	//optimization - steepest ascent w/ line search
	for (int it = 0; it < max_its; it++) {
		//compute the gradient of hitting prob and rate
		hit = computeGradHP(initial, numInteractions, kappaVals, db, particleTypes, 
											Tconst,  ground, targets, gH);
		R = computeGradRate(initial, numInteractions, kappaVals, db, particleTypes, 
												Tconst, targets, gR);

		double obj = hit;
		if (R < c) {
			obj += r*(R-c);
		}


		//apply the rate penalty if needed
		applyRateConstraint(R, c, r, numInteractions, gH, gR);

		//find a step size to give decrease in obj fn
		double step = 0; 
		double delta = lineSearch(initial, numInteractions, kappaVals, db, particleTypes, 
															Tconst, ground, targets, c, r, hit, R, gH, step);

		//take the step
		for (int i = 0; i < numInteractions; i++) {
			kappaVals[i] += step * gH[i];
		}
		checkPositive(numInteractions, kappaVals);

		//output progress
		printf("Iteration: %d, OBJ: %f, Delta: %f\n", it, obj, delta);


		//check for objective function tolerance 
		if (delta < objTol) {
			break;
		}
	}

	std::cout << "Maximum Probability: " << hit << "\n";

	std::cout << kappaVals[0] << ' ' << kappaVals[1] << ' ' << kappaVals[2] << ' ' << "\n";

	double Rf = getRate(initial, kappaVals, db, particleTypes, Tconst, targets);
	printf("The average transition rate with this kappa is %f\n", Rf);


	//free memory
	delete []particleTypes; delete []kappaVals; 
	delete []Tconst; delete []gH; delete []gR; 
}

/****************************************************/
/*********** Some sampling stuff *******************/
/****************************************************/

int findMatrix(int N, double* X, Database* db) {
	//find the adjacency matrix for state X in the database

	//construct the AM
	int* M = new int[N*N]; for (int i = 0; i < N*N; i++) M[i]=0;
	getAdj(X,N,M);
	//refine(N, X, M); //refines X and fills M
	int* dbM = new int[N*N];

	//define the return variable to be -1 if not found
	int state = -1;

	//loop over the database to search
	for (int i = 0; i < db->getNumStates(); i++) {
		extractAM(N, i, dbM, db);
		bool S = checkSame(dbM, M, N);
		if (S) { //found matrix
			state = i;
			break;
		}
	}

	//free memory
	delete []M; delete []dbM;

	return state;
}

void estimateHittingProbability(int N, Database* db, int target) {
	//perform brownian dynamics with file specified interactions to estimate
	//previously determined hitting probabilities

	//set parameters
	int rho = RANGE;
	double beta = BETA;
	int method = 1;
	int num_states = db->getNumStates();
	int num_trajectories = 50;
	int T_cut = 20;
	double DT = 0.01;

	//initialize interaction matrices
	int* types = new int[N];
	int numTypes = bd::readDesignFile(N, types);
	int numInteractions = numTypes*(numTypes+1)/2;
	double* kappa = new double[numInteractions];
	readKappaFile(numInteractions, kappa);
	std::map<std::pair<int,int>, double> kmap; 
	makeKappaMap(numTypes, kappa, kmap);
	int* P = new int[N*N];
	double* E = new double[N*N];
	bd::fillP(N, types, P, E, kmap);

	//set the initial and final state storage
	double* X0 = new double[DIMENSION*N];

	//set potential type
	int pot = POTENTIAL; //0 morse, 1 lj

	//get ground states
	std::vector<int> ground; //vector to hold all ground states
	double* Tconst = new double[num_states*num_states];
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

	//construct a map from state to number of hits
	std::map<int,double> hitSamples;
	for (int i = 0; i < ground.size(); i++) {
		hitSamples.insert(std::pair<int,double>(ground[i],0.0));
	}

	//loop over some number of trajectories to get first hit probs
	int traj = 0;
	while (traj < num_trajectories) {
		setupChain(X0,N);

		double T = 0; 
		while (T < T_cut) {
			//solve the sde, determine new state
			solveSDE(X0, N, DT, rho, beta, E, P, method, pot);
			T += DT;
			int new_state = findMatrix(N, X0, db);

			//check if new state is a ground state, if yes, increment and break
			std::map<int, double>::iterator it = hitSamples.find(new_state); 
			if (it != hitSamples.end()) {
			  it->second ++;
			  traj++;
			  std::cout << "Finished trajectory " << traj << "\n";
			  break;
			}
		}

		if (T > T_cut) {
			std::cout << "Did not find a ground state in the allotted time\n";
		}
	}

	//normalize the empirical probabilities by num_trajectories
	//get total probability for the target
	double targetProb = 0;
	for (auto it = hitSamples.begin(); it != hitSamples.end(); it++) {
		it->second /= num_trajectories;
		double prob = it->second;
		int state = it->first;
		printf("Ground State %d, Probability %f\n", state, prob);

		std::vector<int>::iterator itV = std::find(targets.begin(), targets.end(), state);
		if (itV != targets.end()) {
			targetProb += prob;
		}
	}

	printf("Total Target Probability: %f\n", targetProb);

	//free memory
	delete []types; delete []kappa; delete []P; delete []E; delete []X0;
	delete []Tconst;

}




}