#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "latticeP.h"
#include "../defines.h"

/* Run an mcmc simulation of a lattice protein  */

int main(int argc, char* argv[]) {

	//handle input
	if (argc != 4) {
		fprintf(stderr, "Usage: <Num particles> <db file> <useFile>  %s\n", argv[0]);
		return 1;
	}
	int N = atoi(argv[1]);
	std::string dbFile (argv[2]);
	bool useFile = atoi(argv[3]);

	//get the database here
	lattice::Database* db = lattice::readData(dbFile);

	//test the MCMC functions
	//lattice::runMCMC(N, useFile);

	//build initial database, and update database to get all states
	//lattice::buildPDB(N);
	//lattice::updatePDB(N, db);



	//do mfpt estimation for each state
	int num_states = db->getNumStates();
	for (int i = 0; i < num_states; i++) {
		printf("Beginning estimation for state %d of %d\n", i+1, num_states);
		lattice::estimateMFPT(N, i, db);
	}
	//print the new db
	std::string out = "N" + std::to_string(N) + "mfpt.txt";
	std::ofstream out_str(out);
	out_str << *db;



	//free memory
	delete db;

	return 0;
}