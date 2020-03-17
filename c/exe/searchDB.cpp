#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "bDynamics.h"
#include "database.h"

/* Search a database for a given state by bonds */

int main(int argc, char* argv[]) {

	//handle input
	if (argc < 2) {
		fprintf(stderr, "Usage: <db file> <bond constraint> ... %s\n", argv[0]);
		return 1;
	}
	std::string db_file (argv[1]);
	int bond_cons = atoi(argv[2]);
	int num_opts = argc-3;
	int num_bonds = num_opts / 2;
	std::string* s = new std::string[num_opts];
	for (int i = 0; i < num_opts; i++) {
		s[i] = argv[i+3];
	}

	//get the database here
	bd::Database* db = bd::readData(db_file);
	
	//find states consistent with search requirements
	findState(db, num_bonds, s, bond_cons);


	//free memory 
	delete db; delete []s;

	return 0;
}