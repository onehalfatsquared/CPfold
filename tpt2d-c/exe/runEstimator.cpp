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

	//have a code to filter out un-physical states? todo

	//call the estimator
	bd::estimateMFPT(db->getN(), 620, db);

	//output stuff - for debug
	//std::cout << *db;
	//printf("mfpt = %f\n", (*db)[0].getMFPT());






	//free memory - delete database
	delete db;

	return 0;
}