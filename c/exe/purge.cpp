#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include "bDynamics.h"
#include "database.h"



int main(int argc, char* argv[]) {

	//handle input
	if (argc != 2) {
		fprintf(stderr, "Usage: <Input File> %s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);

	//seed
	srand( time(NULL));

	//get the database here
	bd::Database* db = bd::readData(infile);

	//print messages to user
	printf("Database of states has been read.\n");
	printf("Beginning testing for un-physical states.\n");

	//get the states to purge
	bd::purgeUnphysical(db);

	//purging is done in output
	printf("Writing un-purged states to new file.\n");

	//output stuff - for debug
	//std::cout << *db;

	//output the new states to a file
	std::string out = infile.substr(6,2);
	out = out + "Purge.txt";
	std::ofstream out_str(out);
	out_str << *db; 

	//free memory - delete database
	delete db;

	return 0;
}