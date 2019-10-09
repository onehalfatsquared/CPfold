#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "bDynamics.h"
#include "database.h"

/* Combines two database files given as input and returns a master file */

int main(int argc, char* argv[]) {

	//handle input
	if (argc != 3) {
		fprintf(stderr, "Usage: <Input File 1> <Input File 2>  %s\n", argv[0]);
		return 1;
	}
	std::string infile1 (argv[1]);
	std::string infile2 (argv[2]);

	//get the database here
	bd::Database* db1 = bd::readData(infile1);
	bd::Database* db2 = bd::readData(infile2);

	//compare the two 
	bd::combineMFPTdata(db1, db2);

	//output the new database 
	std::string out = "master.txt";
	std::ofstream out_str(out);
	out_str << *db1; 

	//free memory 
	delete db1; delete db2;

	return 0;
}