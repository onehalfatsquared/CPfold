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

	//lattice::runMCMC(N, useFile);
	//lattice::buildPDB(N);

	lattice::updatePDB(N, db);



	//free memory
	delete db;

	return 0;
}