#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "bDynamics.h"
#include "database.h"

/* Re-writes the filler frequency statistic with states in a database with
   equilibrium probabilities from a seperate sampler in another file */

int main(int argc, char* argv[]) {

	//handle input
	if (argc != 3) {
		fprintf(stderr, "Usage: <Database File> <Eq Probability File>  %s\n", argv[0]);
		return 1;
	}
	std::string infile1 (argv[1]);
	std::string infile2 (argv[2]);

	//get the database here
	bd::Database* db = bd::readData(infile1);

	//update the entries 
	bd::updateFreq(db, infile2);

	//output the new database to a file
	std::string out = infile1;
	out = out + "new.txt";
	std::ofstream out_str(out);
	out_str << *db; 
	std::cout << "Saving new Database to " << out << "\n";

	//free memory 
	delete db; 

	return 0;
}