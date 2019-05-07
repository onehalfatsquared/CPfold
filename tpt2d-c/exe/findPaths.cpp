#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "database.h"
#include "graph.h"



int main(int argc, char* argv[]) {

	//handle input
	if (argc != 3) {
		fprintf(stderr, "Usage: <Input File> <target> %s\n", argv[0]);
		return 1;
	}
	std::string infile (argv[1]);

	int target = atoi(argv[2]);

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

	//print out path with fastest rate
	if (target == -1) {
		QP(g, 1);
	}
	else {
		QP(g, 1, target);
	}

	//construct the subgraph that ends at target
	bd::Graph* sub = bd::targetSubgraph(g, 1, target);

	//print most probable path
	MPP(sub,1);

	//print the graph
	bd::printGraph(sub, 1, 1, 0);

	//free memory - delete database
	delete db; delete g; delete sub;

	return 0;
}