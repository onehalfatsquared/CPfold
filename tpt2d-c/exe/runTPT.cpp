#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include "database.h"
#include "tpt.h"
#include "nauty.h"



int main(int argc, char* argv[]) {

	//handle input
	if (argc < 3 || argc > 4) {
		fprintf(stderr, "Usage: <Input File> <Target State> "
			"Optional <Include Isomorphic>%s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);
	int target = atoi(argv[2]);
	//check for optional input
	bool getIso;
	if (argc == 4) {
		getIso = atoi(argv[3]);
	}
	else {
		getIso = 0;
	}

	//get the database here
	bd::Database* db = bd::readData(infile);

	//specify initial state - worm is 0
	int initial = 0; 

	//run tpt function
	bd::performTPT(db->getN(), initial, target, db, getIso);

	//std::cout << *db;




	//free memory - delete database
	delete db;

	return 0;
}