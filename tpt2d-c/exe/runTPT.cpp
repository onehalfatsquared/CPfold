#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include "database.h"
#include "tpt.h"



int main(int argc, char* argv[]) {

	//handle input
	if (argc != 2) {
		fprintf(stderr, "Usage: <Input File> %s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);

	//get the database here
	bd::Database* db = bd::readData(infile);

	//specify initial and target states. Can be many. 
	int initial; int target;

	//run tpt function
	bd::performTPT(db->getN(), initial, target, db);
	








	//free memory - delete database
	delete db;

	return 0;
}