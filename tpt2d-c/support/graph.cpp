#include "database.h"
#include "graph.h"
#include "pair.h"
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>

namespace bd {

//vertex constructor 
Vertex::Vertex() {
	index = 0; prob = 0; 
}

//vertex deconstructor necceary?//todo
Vertex::~Vertex() {}

//graph constructor
Graph::Graph(int N_) {
	N = N_;
	vertices = new Vertex[N]; 
}

//graph deconstructor
Graph::~Graph() {
	delete []vertices;
}


void updateGraph(int node, Database* db, Graph* g, bool* present, int& count) {
	/*takes a node, makes all its edges, checks if targets exist already,
	creates or updates target nodes*/

	//get the interactions of the current node
	std::vector<Pair> P = (*db)[node].getP();
	double S = (*db)[node].sumP();
	double mfpt = (*db)[node].getMFPT();

	//loop over P - make edges, add nodes
	for (int i = 0; i < P.size(); i++) {
		//make edge
		int target = P[i].index; double prob = P[i].value/S;
		(*g)[node].edges.push_back(Edge(target, 1/mfpt*prob, prob));

		//add/update nodes
		if (!present[target]) {//node does not exist, add
			present[target] = 1; //update present arrays
			Vertex& v = (*g)[target]; //set reference v 
			v.index = target;
			v.prob = v.prob + (*g)[node].prob * prob;
			count++;
		}
		else {//node already exists, just update
			(*g)[target].prob = (*g)[target].prob + (*g)[node].prob * prob;
		}
	}
}


Graph* makeGraph(Database* db) {
	//using the values in db, make a graph from it

	//get the number of particles and states. states = number of nodes in graph
	int N = db->getNumStates(); int particles = db->getN();

	//make an array of node indices already in the graph. 0 ->no, 1->yes.
	bool* present = new bool[N]; for (int i = 0; i < N; i++) present[i] = 0;

	//call graph constructor
	Graph* graph = new Graph(N);

	//set the minimum number of bonds - particles-1
	int min_bond = particles-1; int first = 0;

	//search the database for the worm state and create the vertex
	for (int i = 0; i < N; i++) {
		if ((*db)[i].getBonds() == min_bond) {
			first = i; 
			break;
		}
	}
	int count = 1; //count the number of nodes made - 1 so far
	Vertex& v = (*graph)[first]; //set reference v to first vertex in graph
	v.index = first; //set the index of worm state, just found
	v.prob = 1; //guaranteed starting configuration, prob = 1
	updateGraph(first, db, graph, present, count);


	//loop over states, from top to bottom of the graph, making vertices
	int bond_num = min_bond + 1;
	while (count < N) {
		for (int i = 0; i < N; i++) {
			if ((*db)[i].getBonds() == bond_num) {
				updateGraph(i, db, graph, present, count);
			}
		}
		bond_num++;
	}

	//free memory
	delete []present; 

	//return the graph
	return graph;
}

void endDistribution(Graph* g) {
	//print out the probability distribution for final state

	printf("End Probability Distribution:\n");

	//loop over nodes. CHeck if num edges = 0 -> final state
	for (int i = 0; i < g->getN(); i++) {
		int E = (*g)[i].getNumEdges();
		if (E == 0) {
			printf("State %d, Probability %f\n", i, (*g)[i].getProb());
		}
	}

}

int findMaxProb(Graph* g, int source, int E) {
	//find the edge with the greatest probability

	double maxVal = 0; int maxIndex = 0;
	for (int i = 0; i < E; i++) {
		if ((*g)[source].getEdgeProb(i) > maxVal) {
			maxVal = (*g)[source].getEdgeProb(i);
			maxIndex = (*g)[source].getEdgeTarget(i);
		}
	}

	return maxIndex;
}


void MPP(Graph* g, int source) {
	//determine the most probable path down the graph from given source

	std::vector<int> path; path.push_back(source);
	int E = (*g)[source].getNumEdges();
	while (E > 0) {
		source = findMaxProb(g, source, E);
		path.push_back(source);
		E = (*g)[source].getNumEdges();
	}

	printf("Most Probable Path: ");
	for (int i = 0; i < path.size()-1; i++) {
		printf("%d, ", path[i]);
	}
	printf("%d\n", path[path.size()-1]);


}

void QP(Graph* g, int source) {
	//determine the quickest folding path (largest rate)
	

}










}