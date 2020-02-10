#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "bDynamics.h"
#include "database.h"
#include "../defines.h"

/* Takes in a database with hyperstatic clusters having a 0 equilibrium probability 
   and fills in the values using theory results.  */

int main(int argc, char* argv[]) {

	//handle input
	if (argc != 2) {
		fprintf(stderr, "Usage: <Database File>  %s\n", argv[0]);
		return 1;
	}
	std::string infile1 (argv[1]);

	//get the database here
	bd::Database* db = bd::readData(infile1);

	//set parameters needed for the update
	double beta = BETA;
	double k = KAP;
	double r = RANGE;
	double E = bd::stickyNewton(8.0, r, k, beta);
	printf("When kappa is %f, the well-depth used is %f\n", k, E);

	//update the entries 
	bd::updateHyperstatic(db, beta, E);

	//output the new database to a file
	std::string out = infile1.substr(6,2);
	out = out + "new.txt";
	std::ofstream out_str(out);
	out_str << *db; 

	//free memory 
	delete db; 

	return 0;
}