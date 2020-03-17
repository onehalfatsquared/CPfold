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


bool hasInvalidInteraction(int state, Database* db, int* particleTypes,
													 std::map<std::pair<int,int>,double> kappa) {
	//check if any bonds are of a disallowed type

	double btol = 1;
	int N = db->getN();

	for (int i = 0; i < N; i++) {
		for (int j = i+2; j < N; j++) {
			if ( (*db)[state].isInteracting(i,j,N)) {
				int p1 = particleTypes[i];
				int p2 = particleTypes[j];
				if (kappa[{p1,p2}] < btol) {
					return true;
				}
			}
		}
	}

	return false;

}

void getKappaQuench(double* kappaVals, std::vector<Point>& kappaQuench) {
	//fill a vector with all possible kappa values for a quench
	//this only works for 3 interactions

	Point interaction;

	if (kappaVals[0] == 0) {
		interaction.x = 1; interaction.y = kappaVals[1]; interaction.z = kappaVals[2];
		kappaQuench.push_back(interaction);
	}
	if (kappaVals[1] == 0) {
		interaction.x = kappaVals[0]; interaction.y = 1; interaction.z = kappaVals[2];
		kappaQuench.push_back(interaction);
	}
	if (kappaVals[2] == 0) {
		interaction.x = kappaVals[0]; interaction.y = kappaVals[1]; interaction.z = 1;
		kappaQuench.push_back(interaction);
	}
	if (kappaVals[0] == 0 && kappaVals[1] == 0) {
		interaction.x = 1; interaction.y = 1; interaction.z = kappaVals[2];
		kappaQuench.push_back(interaction);
	}
	if (kappaVals[0] == 0 && kappaVals[2] == 0) {
		interaction.x = 1; interaction.y = kappaVals[1]; interaction.z = 1;
		kappaQuench.push_back(interaction);
	}
	if (kappaVals[1] == 0 && kappaVals[2] == 0) {
		interaction.x = kappaVals[0]; interaction.y = 1; interaction.z = 1;
		kappaQuench.push_back(interaction);
	}
	
}

void getColor(double* kappaVals, std::string& color) {
	//assign a color based on the values in kappaVals

	/* Blue - AA, Red - BB, Yellow - AB, Purple - AA+BB, Orange - AB+BB, Green - AA+BB */

	if (kappaVals[0] > 0 && kappaVals[1] == 0 && kappaVals[2] == 0) {
		color = "blue";
	}
	if (kappaVals[0] == 0 && kappaVals[1] > 0 && kappaVals[2] == 0) {
		color = "yellow";
	}
	if (kappaVals[0] == 0 && kappaVals[1] == 0 && kappaVals[2] > 0) {
		color = "red";
	}
	if (kappaVals[0] > 0 && kappaVals[1] > 0 && kappaVals[2] == 0) {
		color = "green";
	}
	if (kappaVals[0] > 0 && kappaVals[1] == 0 && kappaVals[2] > 0) {
		color = "purple";
	}
	if (kappaVals[0] == 0 && kappaVals[1] > 0 && kappaVals[2] > 0) {
		color = "orange";
	}
	if (kappaVals[0] > 0 && kappaVals[1] > 0 && kappaVals[2] > 0) {
		color = "black";
	}

}

void getPstring(int N, int* particleTypes, std::string& pString) {
	//convert 0 and 1 to A and B

	for (int i = 0; i < N; i++) {
		if (particleTypes[i] == 0) {
			pString.push_back('A');
		}
		else {
			pString.push_back('B');
		}
	}
}

void killTransitions(double* Tnew, int new_states, double* T, int* lumpMap, Database* db, 
										 int* particleTypes, std::map<std::pair<int,int>,double> kappa) {
	//modifies the transition matrix to remove transitions that are no longer 
	//possible under given kappa values

	int N = db->getN();
	int num_states = db->getNumStates();

	bool* toRemove = new bool[num_states];
	//determine the states that need to be removed
	for (int i = 0; i < num_states; i++) {
		if (hasInvalidInteraction(i, db, particleTypes, kappa)) {
			toRemove[i] = true;
		}
		else {
			toRemove[i] = false;
		}

		//std::cout << i << ' ' << toRemove[i] << "\n";
	}

	//set all the removed states transition rates to 0
	for (int i = 0; i < num_states; i++) {
		for (int j = 0; j < num_states; j++) {
			if (toRemove[i] || toRemove[j]) {
				T[toIndex(i,j,num_states)] = 0;
			}

			//condense full matrix to the lumped matrix
			int new_i = lumpMap[i]; int new_j = lumpMap[j];
			Tnew[toIndex(new_i, new_j, new_states)] += T[toIndex(i, j, num_states)];
		}
	}

	//free memory
	delete []toRemove;
}

void killTransitionsQuench(double* Tnew, int new_states, double* T, int* lumpMap, Database* db, 
										 int* particleTypes, std::map<std::pair<int,int>,double> kappa1, 
										 std::map<std::pair<int,int>,double> kappa2, int c1) {
	//modifies the transition matrix to remove transitions that are no longer 
	//possible under given kappa values. This considers a quenched system where the kappa
	//values change at some point

	int N = db->getN();
	int num_states = db->getNumStates();

	//the state that first gets reached is in terms of db2, but we have db1 here. 
	//search the lumpmap for the first entry equal to c1, then get its bonds
	int c1_original;
	for (int i = 0; i < num_states; i++) {
		if (lumpMap[i] == c1) {
			c1_original = i;
		}
	}
	int bond_cut = (*db)[c1_original].getBonds();
	//std::cout << c1 << "\n";
	//std::cout << bond_cut << "\n";

	bool* toRemove = new bool[num_states];
	//determine the states that need to be removed
	for (int i = 0; i < num_states; i++) {
		int bond = (*db)[i].getBonds();
		if (bond <= bond_cut) {
			if (hasInvalidInteraction(i, db, particleTypes, kappa1)) {
				toRemove[i] = true;
			}
			else {
				toRemove[i] = false;
			}
		}
		else{
			if (hasInvalidInteraction(i, db, particleTypes, kappa2)) {
				toRemove[i] = true;
			}
			else {
				toRemove[i] = false;
			}
		}

		//std::cout << i << ' ' << toRemove[i] << "\n";
	}

	//set all the removed states transition rates to 0
	for (int i = 0; i < num_states; i++) {
		for (int j = 0; j < num_states; j++) {
			if (toRemove[i] || toRemove[j]) {
				T[toIndex(i,j,num_states)] = 0;
			}

			//condense full matrix to the lumped matrix
			int new_i = lumpMap[i]; int new_j = lumpMap[j];
			Tnew[toIndex(new_i, new_j, new_states)] += T[toIndex(i, j, num_states)];
		}
	}

	//free memory
	delete []toRemove;
}

void irreversibleGraph(int N, Database* db, int initial, int target, bool useFile) {
	//construct a downhill graph by modifying the interactions via transition matrix 
	// and sending to a graphviz function to be printed

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

	//DO LUMPING - db does not change, just gets assigned a lumpMap
	//upon being output, the new database reflects the lumping
	lumpPerms(db, true);
	std::string out = "temp.txt";
	std::ofstream out_str(out);
	out_str << *db; 
	Database* db2 = readData(out);
	int* lumpMap = new int[num_states];
	for (int i = 0; i < num_states; i++) lumpMap[i] = db->lumpMap[i];

	//make the new arrays
	int new_states = db2->getNumStates();
	double* Fnew = new double[new_states];
	double* eqNew = new double[new_states];
	double* Tnew = new double[new_states*new_states];

	//fill eqNew
	for (int i = 0; i < new_states; i++) {
		eqNew[i] = 0; Fnew[i] = 0;
	}
	for (int old_state = 0; old_state < num_states; old_state++) {
		int new_state = lumpMap[old_state];
		eqNew[new_state] += eq[old_state];
	}

	//compute new free energy
	computeFreeEnergy(new_states, eqNew, Fnew);

	//zero out new transition matrix
	for (int i = 0; i < new_states*new_states; i++) {
		Tnew[i] = 0;
	}
	//make the Tnew matrix - zero out any rates to inaccesible state
	killTransitions(Tnew, new_states, T, lumpMap, db, particleTypes, kappa);

	//make a graph structure of the database
	int node;
	double prob;
	Graph* g = makeGraphTM(Tnew, initial, N, new_states, node, prob);

	//print out graphviz file with tpt data
	printGraphPF(g, initial, Fnew, 1, 1, 1);

	//free the memory
	delete Z;
	delete []T; delete []eq; delete []F; 
	delete []Tnew; delete []eqNew; delete []Fnew;
	delete []particleTypes; delete []kappaVals;
	delete []lumpMap;
	delete db2; delete g;
}

bool generateGraph(int N, Database* db, int initial, int numTypes, int numInteractions, 
									 std::vector<int> ground, int* particleTypes, Point interaction, 
									 double* Tconst, int& max_node, bool final) {
	/*using the given particle types and interaction types, creates the 
	modified transition matrix, makes a graph, and checks if it ends in a 
	node with probability 1. */

	//if final is true,only return success if a ground state is reached

	//get number of states
	int num_states = db->getNumStates();

	//store transition matrix in new array
	double* T = new double[num_states*num_states];
	for (int i = 0; i < num_states*num_states; i++) {
		T[i] = Tconst[i];
	}

	//set the kappa values
	double* kappaVals = new double[numInteractions];
	kappaVals[0] = interaction.x; 
	kappaVals[1] = interaction.y; 
	kappaVals[2] = interaction.z; 

	//init the equilibrium measure
	double* eq = new double[num_states];  

	//declare the kappa map
	std::map<std::pair<int,int>,double> kappa;

	//compute the eq probability at the current kappa vals
	//make the map
	makeKappaMap(numTypes, kappaVals, kappa);
	//do rewieght
	reweight(N, num_states, db, particleTypes, eq, kappa);

	//DO LUMPING - db does not change, just gets assigned a lumpMap
	//upon being output, the new database reflects the lumping
	lumpPerms(db, true);
	std::string out = "temp.txt";
	std::ofstream out_str(out);
	out_str << *db; 
	Database* db2 = readData(out);
	int* lumpMap = new int[num_states];
	for (int i = 0; i < num_states; i++) lumpMap[i] = db->lumpMap[i];

	//make the new arrays
	int new_states = db2->getNumStates();
	double* Tlump = new double[new_states*new_states];

	//zero out new transition matrix
	for (int i = 0; i < new_states*new_states; i++) {
		Tlump[i] = 0;
	}
	//make the Tnew matrix - zero out any rates to inaccesible state
	killTransitions(Tlump, new_states, T, lumpMap, db, particleTypes, kappa);

	//make a graph structure of the database
	int node;
	double prob;
	Graph* g = makeGraphTM(Tlump, initial, N, new_states, node, prob);

	//determine if this graph should be output -  if prob is close to 1
	bool success = false; 
	max_node = 0;
	int bonds = (*db2)[node].getBonds();
	if (prob > 0.98 && !final && bonds >= N+1) {
		success = true;
		max_node = node;
	}
	if (prob > 0.98 && final) {
		if (bonds >= 2*N-3) {
			success = true;
			max_node = node;
		}
	}

	//free memory
	delete []T; delete []eq; delete []lumpMap; delete []Tlump;
	delete db2; delete g; delete []kappaVals;

	//return the bool
	return success;
}

void generateAllGraphs(int N, Database* db, int initial) {
	//generate all graphs that show the irreversible folding of a linear chain
	//Considers all unique permutations of particle types
	//Considers all unique combinations of on/off interactions
	//hardcoded to 2 types of particles and 3 interactions

	//if a graph ends in a node of 100% probability, save those settings in a file
	//file order - cluster id - particle types - kappa values

	//get database info
	int num_states = db->getNumStates(); 

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
	std::vector<int> ground; //vector to hold all ground states - also target states here
	createTransitionMatrix(Tconst, num_states, db, ground);

	//get all the permutations of particle types
	std::deque<std::string> particle_perms; 
	distinctPerms(N, particle_perms);
	int num_particle_perms = particle_perms.size();

	//get all the permutations of interactions
	std::vector<Point> interaction_perms;
	int num_interaction_perms = 6;
	interaction_perms.push_back(Point(10, 0, 0));
	interaction_perms.push_back(Point(0, 10, 0));
	interaction_perms.push_back(Point(0, 0, 10));
	interaction_perms.push_back(Point(10, 10, 0));
	interaction_perms.push_back(Point(10, 0, 10));
	interaction_perms.push_back(Point(0, 10, 10));

	//open a file to write the data to
	std::ofstream ofile;
	ofile.open("possibleQuenching.txt");
	ofile << N << "\n";

	//loop over each permutation of particles
	for (int pperm = 0; pperm < num_particle_perms; pperm++) {
		//set the current permutation in particleTypes
		setTypes(N, particleTypes, particle_perms, pperm);
		printf("Testing particle permutation %d of %d\n", pperm+1, num_particle_perms);

		//loop over the interaction perms
		for (int iperm = 0; iperm < num_interaction_perms; iperm++) {
			//call a function to construct the transition matrix, make the
			//graph, and check if it should be output

			int max_node = initial;
			Point interaction = interaction_perms[iperm];
			bool success = generateGraph(N, db, initial, numTypes, numInteractions, ground, particleTypes, 
																	 interaction, Tconst, max_node, false);

			if (success) {
				printf("Cluster %d is a guaranteed endpoint\n", max_node);
				//output the cluster id
				ofile << max_node << ' ';
				//output the particle type distribution
				for (int i = 0; i < N; i++) {
					ofile << particleTypes[i] << ' ';
				}
				//output the interaction matrix
				ofile << interaction.x << ' ';
				ofile << interaction.y << ' ';
				ofile << interaction.z << ' ';
				//make new line
				ofile << "\n";
			}
		}
	}

	//close the outfile
	ofile.close();

	//free the memory
	delete []particleTypes; delete []Tconst;
}


void testQuench(int N, Database* db, int initial, int numTypes, int numInteractions, 
								std::vector<int> ground, int* particleTypes, std::ofstream& ofile, 
								double* Tconst, double* kappaVals) {
	//determines if a quenching scheme to a definite ground state exists

	//must construct all quenched systems from the given kappaVals
	std::vector<Point> kappaQuench;
	getKappaQuench(kappaVals, kappaQuench);

	//loop over all possible quenches, do the test
	for (int i = 0; i < kappaQuench.size(); i++) {
		Point interaction = kappaQuench[i];

		//generate graph, check if it ends in ground state w/ probability 1
		int max_node = 0;
		bool success = generateGraph(N, db, initial, numTypes, numInteractions, ground, 
															   particleTypes, interaction, Tconst, max_node, true);

		//if yes, print out to new file
		if (success) {
			printf("Cluster %d is a guaranteed endpoint\n", max_node);
			//output the particle type distribution
			for (int i = 0; i < N; i++) {
				ofile << particleTypes[i] << ' ';
			}
			//output the initial cluster
			ofile << initial << ' ';
			//output the original interactions
			ofile << kappaVals[0] << ' ';
			ofile << kappaVals[1] << ' ';
			ofile << kappaVals[2] << ' ';
			//output the final cluster
			ofile << max_node << ' ';
			//output the final interaction matrix
			ofile << interaction.x << ' ';
			ofile << interaction.y << ' ';
			ofile << interaction.z << ' ';
			//make new line
			ofile << "\n";
		}
	}

	
}

void findQuenches(int N, Database* db) {
	//use the possibleQuenches file to try and reach a ground state

	//get database info
	int num_states = db->getNumStates(); 

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
	std::vector<int> ground; //vector to hold all ground states - also target states here
	createTransitionMatrix(Tconst, num_states, db, ground);

	//open a file to write the new data to
	std::ofstream ofile;
	ofile.open("quenches.txt");
	ofile << N << "\n";

	//open the input file with all the starting points
	std::string filename = "possibleQuenching.txt";
	std::ifstream in_str(filename);

	//read in the input and loop over all possibilities
	int temp; in_str >> temp; //get rid of N at beginning
	int clusterID;

	//loop stores cluster id
	while (in_str >> clusterID) {
		//store particle types
		for (int i = 0; i < N; i++) {
			in_str >> particleTypes[i];
		}

		//store kappa vals
		for (int i = 0; i < numInteractions; i++) {
			in_str >> kappaVals[i];
		}

		//checks for any quench of kappavals that give prob 1 ground state
		testQuench(N, db, clusterID, numTypes, numInteractions, 
								ground, particleTypes, ofile, Tconst, kappaVals);

	}


	//close files
	ofile.close();
	in_str.close();

	//free memory
	delete []particleTypes; delete []kappaVals; delete []Tconst;
}

void generateQuenchedGraph(int N, Database* db, int initial, int numTypes, int numInteractions, 
									 std::vector<int> ground, int* particleTypes, double* Tconst, 
									 double* kappaVals1, double* kappaVals2, int c1, int count) {
	/*generate a graph for the quenched system using the specified parameters */

	//get number of states
	int num_states = db->getNumStates();

	//store transition matrix in new array
	double* T = new double[num_states*num_states];
	for (int i = 0; i < num_states*num_states; i++) {
		T[i] = Tconst[i];
	}

	//init the equilibrium measure
	double* eq = new double[num_states];  

	//declare the kappa map
	std::map<std::pair<int,int>,double> kappa1;
	std::map<std::pair<int,int>,double> kappa2;

	//compute the eq probability at the current kappa vals
	//make the map
	makeKappaMap(numTypes, kappaVals1, kappa1);
	makeKappaMap(numTypes, kappaVals2, kappa2);
	//do rewieght -  is this necessary?
	reweight(N, num_states, db, particleTypes, eq, kappa1);

	//DO LUMPING - db does not change, just gets assigned a lumpMap
	//upon being output, the new database reflects the lumping
	lumpPerms(db, true);
	std::string out = "temp.txt";
	std::ofstream out_str(out);
	out_str << *db; 
	Database* db2 = readData(out);
	int* lumpMap = new int[num_states];
	for (int i = 0; i < num_states; i++) lumpMap[i] = db->lumpMap[i];

	//make the new arrays
	int new_states = db2->getNumStates();
	double* Tlump = new double[new_states*new_states];

	//zero out new transition matrix
	for (int i = 0; i < new_states*new_states; i++) {
		Tlump[i] = 0;
	}
	//make the Tnew matrix - zero out any rates to inaccesible state
	killTransitionsQuench(Tlump, new_states, T, lumpMap, db, particleTypes, kappa1, kappa2, c1);

	//make a graph structure of the database
	int node;
	double prob;
	Graph* g = makeGraphTM(Tlump, initial, N, new_states, node, prob);

	//print out graphviz file 
	int layer_swap = (*db2)[c1].getBonds() - (N-1);
	std::string color1; getColor(kappaVals1, color1);
	std::string color2; getColor(kappaVals2, color2);
	//std::cout << color1 << ' ' << color2 << "\n";
	//convert 0 and 1 types to A and B
	std::string pString;
	getPstring(N, particleTypes, pString);
	printGraphQuenched(g, initial, count, layer_swap, color1, color2, pString, 1, 1, 1);
	

	//free memory
	delete []T; delete []eq; delete []lumpMap; delete []Tlump;
	delete db2; delete g;
}

void graphQuenches(int N, Database* db, int initial) {
	//graph the quenches found. needs to seperate the transition matrices to get the ful path
	//use different colors and the combos to represent bond types to make path possible
	/* Blue - AA, Red - BB, Yellow - AB, Purple - AA+BB, Orange - AB+BB, Green - AA+BB */

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes = 2;
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals1 = new double[numInteractions];
	double* kappaVals2 = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states - also target states here
	createTransitionMatrix(Tconst, num_states, db, ground);

	//open the input file with all the starting points
	std::string filename = "quenches.txt";
	std::ifstream in_str(filename);

	//read in the input and loop over all possibilities
	int temp; in_str >> temp; //get rid of N at beginning
	int clusterMID; int clusterEND;

	//loop stores cluster id
	int count = 0;
	while (in_str >> particleTypes[0]) {
		//store particle types
		for (int i = 1; i < N; i++) {
			in_str >> particleTypes[i];
		}
		//store first cluster hit
		in_str >> clusterMID;

		//store initial kappa vals
		for (int i = 0; i < numInteractions; i++) {
			in_str >> kappaVals1[i];
		}

		//store final cluster hit
		in_str >> clusterEND;

		//store final kappa vals
		for (int i = 0; i < numInteractions; i++) {
			in_str >> kappaVals2[i];
		}

		//debug stuff
		/*
		std::cout << particleTypes[0] << particleTypes[1] << particleTypes[2] << 
		particleTypes[3] << particleTypes[4] << particleTypes[5] << "\n";
		std::cout << clusterMID << "\n";
		std::cout << kappaVals1[0] << kappaVals1[1] <<  kappaVals1[2] << "\n";
		std::cout << clusterEND << "\n";
		std::cout << kappaVals2[0] << kappaVals2[1] <<  kappaVals2[2] << "\n";
		*/

		count++;
		generateQuenchedGraph(N, db, initial, numTypes, numInteractions, ground, 
													particleTypes, Tconst, kappaVals1, kappaVals2, clusterMID, count);
		
	}


	//close files
	in_str.close();

	//free memory
	delete []particleTypes; delete []kappaVals1; delete []kappaVals2; delete []Tconst;
}
















}