#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "bDynamics.h"
#include "database.h"
#include "hydro.h"
#include "../defines.h"

/* Takes in cluster data from hydrodynamics sim and constructs the transition matrix.  */
// Read in an empty database (in purge folder) and write a new db file from this data

int main(int argc, char* argv[]) {

	//handle input
	if (argc != 5) {
		fprintf(stderr, "Usage: <Hydrodynamics data> <empty db file> <dt> <n_save>  %s\n", argv[0]);
		return 1;
	}

	//read input and set parameters
	std::string infile (argv[1]);
	std::string infile2 (argv[2]);
	double dt = atof(argv[3]);
	int n_save = atoi(argv[4]);
	int maxT = 2000;

	//construct empty database from file with no mfpt data
	bd::Database* db = bd::readData(infile2);
	std::cout << "Empty database has been created. \n";

	//get num of particles from db
	int N = db->getN();


	//fill data structure with hydrodynamics data
	std::cout << "Reading in Hydrodynamics data. \n";
	bd::HCC* hc = bd::extractData(infile, N, maxT);
	std::cout << "Hydro data has been stored. \n";

	//bd::testExtract(hc);
	bd::determineTransitions(hc, db);

	//get the transition data - in db

	//output the database to file





	

	delete hc; delete db;
	return 0;
}