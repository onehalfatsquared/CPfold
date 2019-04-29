#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "bDynamics.h"
#include "database.h"

//MAKE SURE TO RUN THIS ON A PURGED DATA SET EX input/N7Purge.txt

int main(int argc, char* argv[]) {

	//handle input
	if (argc != 2) {
		fprintf(stderr, "Usage: <Input File> %s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);

	//get the database here
	bd::Database* db = bd::readData(infile);

	//set database variables
	int N = db->getN(); int num_states = db->getNumStates();

	//tell user information
	printf("Database of states has been read.\n");
	printf("Mean first passage time estimator beginning.\n");


	//call estimator over every state, i. if i has 11 bonds (N=7), do nothing.
	for (int i = 0; i < num_states; i++) {
		if (((*db)[i].getBonds() >= 11 && N == 7) || ((*db)[i].getBonds() >= 9 && N == 6)) {//these states are rigid
			//do nothing
		}
		else{
			//call the estimator
			bd::estimateMFPT(N, i, db);
		}
	}

	//output stuff - for debug
	//std::cout << *db;

	//output the mfpt results to file
	std::string out = infile.substr(6,2);
	out = out + "mfpt.txt";
	std::ofstream out_str(out);
	out_str << *db; 

	//free memory - delete database
	delete db;

	return 0;
}