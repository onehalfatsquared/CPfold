#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include "database.h"
#include "tpt.h"
#include "nauty.h"



int main(int argc, char* argv[]) {

	//handle input
	if (argc < 4 || argc > 5) {
		fprintf(stderr, "Usage: <Input File> <Initial State> <Target State> "
			"Optional <Include Isomorphic>%s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);
	int initial = atoi(argv[2]);
	int target = atoi(argv[3]);
	//check for optional input
	bool getIso;
	if (argc == 5) {
		getIso = atoi(argv[4]);
	}
	else {
		getIso = 0;
	}

	//get the database here
	bd::Database* db = bd::readData(infile);

	//run tpt function
	bd::performTPT(initial, target, db, getIso);

	//get hitting probs
	bd::getHittingProbabilityGS(initial, db);

	bd::makeProbCompareGraph(initial, db);

	//std::cout << *db;




	//free memory - delete database
	delete db;

	return 0;
}