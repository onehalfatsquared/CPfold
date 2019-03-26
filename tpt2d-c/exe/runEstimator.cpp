#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include "bDynamics.h"
#include "database.h"



int main(int argc, char* argv[]) {

	//handle input
	if (argc != 2) {
		fprintf(stderr, "Usage: <Input File> %s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);

	//get the database here
	bd::Database* db = bd::readData(infile);

	//call the estimator
	bd::estimateMFPT(db->getN(), 111, db);

	//output stuff
	//std::cout << *db;






	//free memory - delete database
	delete db;

	return 0;
}