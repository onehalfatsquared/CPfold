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

void lumpArrays(int old_states, int new_states, int* lumpMap, Database* db, double* F, 
								double* Fnew, double* flux, double* fluxNew, double* r, double* rNew, 
								double* eq, double* eqNew, double* q, double* qNew, double* Tnew, 
								int initial, int target) {
	//uses lumpMap to combine the additive quantities measured for the full system

	//first do all the additive quantities - r, eq, q
	//begin by zeroing out all the arrays
	for (int i = 0; i < new_states; i++) {
		rNew[i] = eqNew[i] = qNew[i] = 0;
	}
	//next add all the old values into the lumped versions
	for (int old_state = 0; old_state < old_states; old_state++) {
		int new_state = lumpMap[old_state];
		rNew[new_state] += r[old_state];
		eqNew[new_state] += eq[old_state];
		//qNew[new_state] += q[old_state];
	}

	//compute free energy
	computeFreeEnergy(new_states, eqNew, Fnew);

	//next do the matrix quantities
	for (int i = 0; i < new_states*new_states; i++) {
		Tnew[i] = fluxNew[i] = 0;
	}
	//make the Tnew matrix
	std::vector<int> groundStates;
	createTransitionMatrix(Tnew, new_states, db, groundStates);
	//fill in transposed entries such that Tnew satisfies detailed balance
	satisfyDB(Tnew, new_states, db, eqNew);
	//fill in diagonal with negative sum of row entries
	fillDiag(Tnew, new_states);

	//compute comittor
	std::vector<int> targets; targets.push_back(lumpMap[target]);
	computeCommittor(qNew, Tnew, new_states, lumpMap[initial], targets);

	//make the new flux matrix
	computeFlux(new_states, qNew, Tnew, eqNew, fluxNew);
}

void performTPT(int N, Database* db, int initial, int target, bool useFile) {
	//perform tpt calculations from initial to target states

	//Note: PERMUTATIONS MUST BE LUMPED FOR GRAPHICAL OUTPUT
	//DO the reweighting first, then lump everything according to lumpMap

	//parameters
	int num_states = db->getNumStates();
	std::vector<int> groundStates;

	//init and compute configurational partition function and free energy
	double* Z = new double[num_states]; double* F = new double[num_states];
	for (int i = 0; i < num_states; i++) Z[i] = F[i] = 0;
	computePartitionFn(num_states, db, Z); 
	computeFreeEnergy(num_states, Z, F);

	//initialize the transition rate matrix
	double* T = new double[num_states*num_states];
	for (int i = 0; i < num_states*num_states; i++) {
		T[i] =  0;
	}

	//get bonds -> bonds+1 entries from mfpt estimates
	createTransitionMatrix(T, num_states, db, groundStates);

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
	readKappaFile(numInteractions, kappaVals);

	//init the equilibrium measure
	double* eq = new double[num_states];  

	//declare the kappa map
	std::map<std::pair<int,int>,double> kappa;

	//compute the eq probability at the current kappa vals
	//make the map
	makeKappaMap(numTypes, kappaVals, kappa);
	//do rewieght
	reweight(N, num_states, db, particleTypes, eq, kappa);

	//find all target states consistent with input target
	std::vector<int> targets;
	findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

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

	//init and compute the probability fluxes
	double* flux = new double[num_states*num_states]; 
	for (int i = 0; i < num_states*num_states; i++) flux[i] = 0;
	computeFlux(num_states, q, T, eq, flux);

	//compute rates
	double R = computeTransitionRateTPT(num_states, q, T, eq);
	std::cout << "TPT rate " << R << "\n";
	double* r = new double[num_states];
	computeMFPTs(num_states, T, targets, r);
	//go from mfpt to rates
	for (int i = 0; i < num_states; i++) r[i] = 1.0 / r[i];

	//DO LUMPING - db does not change, just gets assigned a lumpMap
	//upon being output, the new database reflects the lumping
	lumpPerms(db);
	std::string out = "temp.txt";
	std::ofstream out_str(out);
	out_str << *db; 
	Database* db2 = readData(out);
	int* lumpMap = new int[num_states];
	for (int i = 0; i < num_states; i++) lumpMap[i] = db->lumpMap[i];

	//make the new arrays
	int new_states = db2->getNumStates();
	double* Fnew = new double[new_states];
	double* rNew = new double[new_states];
	double* qNew = new double[new_states];
	double* eqNew = new double[new_states];
	double* fluxNew = new double[new_states*new_states];
	double* Tnew = new double[new_states*new_states];

	//call the function to set the lumped variables
	lumpArrays(num_states, new_states, lumpMap, db2, F, Fnew, flux, fluxNew, r, rNew, 
						 eq, eqNew, q, qNew, Tnew, initial, target);
	std::cout << "Average transition rate is " << rNew[initial] << "\n";
	

	//make a graph structure of the database
	Graph* g = makeGraph(db2);
	//std::cout << "graph made";

	//print out graphviz file with tpt data
	printGraphRev(g, initial, Fnew, fluxNew, 1, 1, 0);

	//free the memory
	delete Z;
	delete []T; delete []q; delete []eq; delete []flux; delete []F; delete []r;
	delete []Tnew; delete []qNew; delete []eqNew; delete []fluxNew; 
	delete []Fnew; delete []rNew;
	delete []particleTypes; delete []kappaVals;
	delete []lumpMap;
	delete db2; delete g;

}































}