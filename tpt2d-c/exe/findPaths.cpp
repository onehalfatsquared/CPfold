#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "database.h"
#include "graph.h"



int main(int argc, char* argv[]) {

	//handle input
	if (argc != 3) {
		fprintf(stderr, "Usage: <Input File> <path>%s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);

	int path = atoi(argv[2]);

	//get the database here
	bd::Database* db = bd::readData(infile);
	printf("Database of states has been read.\n");
	
	//make the graph
	bd::Graph* g = bd::makeGraph(db);
	printf("Database has been stored as a graph.\n");

	//print out the ending distribution
	endDistribution(g);

	//print out most probable path
	MPP(g, 1);

	//free memory - delete database
	delete db; delete g;

	return 0;
}