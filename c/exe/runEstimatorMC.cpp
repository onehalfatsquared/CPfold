#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "../defines.h"
#include "bDynamics.h"
#include "database.h"
#include "sampling.h"

//MAKE SURE TO RUN THIS ON A PURGED DATA SET EX input/N7Purge.txt
// Run Type: 0 for full run
//					 1 for chain estimate, with resets
// Chain State: index of the linear chain in the DB. 1 for 6, 0 for 7.

int main(int argc, char* argv[]) {

	//handle input
	if (argc < 3) {
		fprintf(stderr, "Usage: <Input File> <Run Type> <Chain State> %s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);
	int rType = atoi(argv[2]); 
	int source;
	if (rType == 1) {
		if (argc != 4) {
			fprintf(stderr, "Usage: <Input File> <Run Type> <Chain State> %s\n", argv[0]);
			return 1;
		}
		source = atoi(argv[3]); 
	} 

	//get the database here
	bd::Database* db = bd::readData(infile);

	//set database variables
	int N = db->getN(); int num_states = db->getNumStates();

	//initialize the rng
	RandomNo* rngee = new RandomNo();

	//tell user information - estimator beginning
	printf("Database of states has been read.\n");
	printf("Mean first passage time estimator beginning.\n");


	if (rType == 0) {
		//do a full run over all states in DB

		//call estimator over every state, i. if i has 11 bonds (N=7), do nothing.
		for (int i = 0; i < num_states; i++) {
#if (DIMENSION == 2)
			if (((*db)[i].getBonds() >= 11 && N == 7) || ((*db)[i].getBonds() >= 9 && N == 6)) {//these states are rigid
				//do nothing
			}
#elif (DIMENSION == 3) 
			if ((*db)[i].getBonds() >= 12 && N == 6) {//these states are rigid
				//do nothing
			}
#endif
			else{
				//call the estimator
				mcm::estimateMFPTreflect(N, i, db, rngee);
				//mcm::estimateMFPTreset(N, i, db, rngee);
			}
		}

		//output stuff - for debug
		//std::cout << *db;

		//output the mfpt results to file
		std::string out = infile.substr(6,2);
		out = out + "mfptMCTEST2.txt";
		std::ofstream out_str(out);
		out_str << *db; 

	}

	else if (rType == 1) {
		//just perform test on linear chain with resets

		/*
		bd::estimateChain(N, source, db);

		//output the mfpt results to file
		std::string out = infile.substr(6,2);
		out = out + "mfptChainMC.txt";
		std::ofstream out_str(out);
		out_str << *db; 
		*/
	}

	//free memory 
	delete db; delete rngee;

	return 0;
}