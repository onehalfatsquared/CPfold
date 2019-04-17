#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "database.h"
#include "graph.h"



int main(int argc, char* argv[]) {

	//handle input
	if (argc != 4) {
		fprintf(stderr, "Usage: <Input File> <draw> <reduce> %s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);

	int draw = atoi(argv[2]);
	int reduce = atoi(argv[3]);

	//get the database here
	bd::Database* db = bd::readData(infile);

	//print messages to user
	printf("Database of states has been read.\n");
	printf("Constructing GraphViz Code.\n");

	//make graphviz code
	bd::makeGraphViz(db, draw, reduce);

	//free memory - delete database
	delete db;

	return 0;
}