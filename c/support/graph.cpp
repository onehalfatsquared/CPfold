#include "database.h"
#include "graph.h"
#include "pair.h"
#include <map>
#include <vector>
#include <deque>
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
	end = 0;
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
	double tol = 1e-5;

	//loop over P - make edges, add nodes
	for (int i = 0; i < P.size(); i++) {
		//make edge
		int target = P[i].index; double prob = P[i].value/S;
		if (1/mfpt*prob > tol) {
			(*g)[node].edges.push_back(Edge(target, 1/mfpt*prob, prob));
			(*g)[target].back_edges.push_back(node);

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

	//set counters for number of end states
  int countUpdates = 1;

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
				countUpdates += 1;
				//printf("%d done, count = %d\n", i, count);
			}
		}
		bond_num++; 
	}

	//set number of end states
	graph->end = N - countUpdates;
	printf("End states: %d.\n", graph->end);

	//free memory
	delete []present; 

	//return the graph
	return graph;
}

Graph* targetSubgraph(Graph* g, int source, int target) {
	//constructs the subgraph of g corresponding of all states that can reach target

	//first get a list of nodes to include by performing bfs on reverse graph
	int ns = g->getN();
	bool* inSub = new bool[ns]; for (int i = 0; i < ns; i++) inSub[i]=0;
	inSub[target] = 1;
	std::deque<int> toVisit; toVisit.push_back(target);
	while (!toVisit.empty()) {
		int node = toVisit[0];
		int E = (*g)[node].getNumBackEdges();
		for (int edge = 0; edge < E; edge++) {
			int s = (*g)[node].getEdgeSource(edge);
			if (inSub[s] == 0) {
				inSub[s] = 1;
				toVisit.push_back(s);
			}
		}
		toVisit.pop_front();
	}

	//make a new graph that only includes states in inSub

	//present array to check if a state has been added yet
	bool* present = new bool[ns]; for (int i = 0; i < ns; i++) present[i] = 0;

	//call graph constructor
	Graph* graph = new Graph(ns);

	//make the source node
	Vertex& v = (*graph)[source];
	v.index = source;
	v.prob = 1;
	present[source] = 1;

	//loop over source edges. get new normalizer, make edges. add new nodes.
	int E = (*g)[source].getNumEdges();
	double norm = 0;
	for (int i = 0; i < E; i++) {
		int T = (*g)[source].getEdgeTarget(i);
		if (inSub[T]) {
			norm += (*g)[source].getEdgeRate(i);
		}
	}
	for (int i = 0; i < E; i++) {
		int T = (*g)[source].getEdgeTarget(i);
		if (inSub[T]) {
			double rate = (*g)[source].getEdgeRate(i);
			double prob = rate / norm;
			(*graph)[source].edges.push_back(Edge(T, rate, prob));
			(*graph)[T].back_edges.push_back(source);

			if (!present[T]) {//node does not exist, add
				present[T] = 1; //update present arrays
				toVisit.push_back(T);
				Vertex& v = (*graph)[T]; //set reference v 
				v.index = T;
				v.prob = v.prob + (*graph)[source].prob * prob;
				}
			else {//node already exists, just update
				(*graph)[T].prob = (*graph)[T].prob + (*graph)[source].prob * prob;
			}
		}
	}

	//loop over vertices, top down in toVisit, check if should be present, perform update
	while(!toVisit.empty()) {
		int node = toVisit[0];
		int E = (*g)[node].getNumEdges();
		double norm = 0;
		for (int i = 0; i < E; i++) {
			int T = (*g)[node].getEdgeTarget(i);
			if (inSub[T]) {
				norm += (*g)[node].getEdgeRate(i);
			}
		}
		for (int i = 0; i < E; i++) {
			int T = (*g)[node].getEdgeTarget(i);
			if (inSub[T]) {
				double rate = (*g)[node].getEdgeRate(i);
				double prob = rate / norm;
				(*graph)[node].edges.push_back(Edge(T, rate, prob));
				(*graph)[T].back_edges.push_back(node);

				if (!present[T]) {//node does not exist, add
					present[T] = 1; //update present arrays
					toVisit.push_back(T);
					Vertex& v = (*graph)[T]; //set reference v 
					v.index = T;
					v.prob = v.prob + (*graph)[node].prob * prob;
					}
				else {//node already exists, just update
					(*graph)[T].prob = (*graph)[T].prob + (*graph)[node].prob * prob;
				}
			}
		}
		toVisit.pop_front();
	}


	//graph is made. free memory and return graph.
	delete []inSub; delete []present;
	return graph;
}

Graph* sourceSubgraph(Graph* g, int source) {
	//construct the subgraph of g that starts at node source

	//create new graph structure
	int ns = g->getN();
	Graph* graph = new Graph(ns);
	
	//present array to check if a state has been added yet
	bool* present = new bool[ns]; for (int i = 0; i < ns; i++) present[i] = 0;

	//make the source node
	Vertex& v = (*graph)[source];
	v.index = source;
	v.prob = 1;
	present[source] = 1;

	std::deque<int> toVisit; 

	//loop over source edges. get new normalizer, make edges. add new nodes.
	int E = (*g)[source].getNumEdges();
	double norm = 0;
	for (int i = 0; i < E; i++) {
		int T = (*g)[source].getEdgeTarget(i);
		norm += (*g)[source].getEdgeRate(i);
	}
	for (int i = 0; i < E; i++) {
		int T = (*g)[source].getEdgeTarget(i);
		double rate = (*g)[source].getEdgeRate(i);
		double prob = rate / norm;
		(*graph)[source].edges.push_back(Edge(T, rate, prob));
		(*graph)[T].back_edges.push_back(source);

		if (!present[T]) {//node does not exist, add
			present[T] = 1; //update present arrays
			toVisit.push_back(T);
			Vertex& v = (*graph)[T]; //set reference v 
			v.index = T;
			v.prob = v.prob + (*graph)[source].prob * prob;
			}
		else {//node already exists, just update
			(*graph)[T].prob = (*graph)[T].prob + (*graph)[source].prob * prob;
		}
	}

	//loop over vertices, top down in toVisit, check if should be present, perform update
	while(!toVisit.empty()) {
		int node = toVisit[0];
		int E = (*g)[node].getNumEdges();
		double norm = 0;
		for (int i = 0; i < E; i++) {
			int T = (*g)[node].getEdgeTarget(i);
			norm += (*g)[node].getEdgeRate(i);
		}
		for (int i = 0; i < E; i++) {
			int T = (*g)[node].getEdgeTarget(i);
			double rate = (*g)[node].getEdgeRate(i);
			double prob = rate / norm;
			(*graph)[node].edges.push_back(Edge(T, rate, prob));
			(*graph)[T].back_edges.push_back(node);

			if (!present[T]) {//node does not exist, add
				present[T] = 1; //update present arrays
				toVisit.push_back(T);
				Vertex& v = (*graph)[T]; //set reference v 
				v.index = T;
				v.prob = v.prob + (*graph)[node].prob * prob;
				}
			else {//node already exists, just update
				(*graph)[T].prob = (*graph)[T].prob + (*graph)[node].prob * prob;
			}
		}
		toVisit.pop_front();
	}


	//graph is made. free memory and return graph.
	delete []present;
	return graph;
}

void endDistribution(Graph* g) {
	//print out the probability distribution for final state

	printf("End Probability Distribution:\n");

	//loop over nodes. CHeck if num edges = 0 -> final state
	for (int i = 0; i < g->getN(); i++) {
		int E = (*g)[i].getNumEdges();
		int B = (*g)[i].getNumBackEdges();
		if (E == 0 && B != 0) {
			printf("State %d, Probability %f\n", i, (*g)[i].getProb());
		}
	}
}

void tOrder(Graph* g, int* T, int source) {
	//construct a topological ordering of the graph g

	std::deque<int> queue; std::vector<int> done;
	T[0] = source; queue.push_back(source);
	int count = 1;
	while (!queue.empty()) {
		int node = queue[0];
		int edges = (*g)[node].getNumEdges();
		for (int i = 0; i < edges; i++) {
			int new_node = (*g)[node].getEdgeTarget(i);
			if (std::find(done.begin(), done.end(), new_node) == done.end()) {
				T[count] = new_node; queue.push_back(new_node); 
				count++; done.push_back(new_node);
			}
		}
		queue.pop_front();
	}
}

void findMaxProb(Graph* g, int& source, int E, std::vector<int>& path, std::vector<double>& prob) {
	//find the edge with the greatest probability

	double maxVal = 0; int maxIndex = 0;
	for (int i = 0; i < E; i++) {
		if ((*g)[source].getEdgeProb(i) > maxVal) {
			maxVal = (*g)[source].getEdgeProb(i);
			maxIndex = (*g)[source].getEdgeTarget(i);
		}
	}
	path.push_back(maxIndex); prob.push_back(maxVal);
	source = maxIndex;
}



void fillRates(Graph* g, std::vector<int> path, std::vector<double>& rates) {
	//determine the rate for each step of a given path

	for (int i = 0; i < path.size()-1; i++) {
		int source = path[i]; int target = path[i+1];
		int E = (*g)[source].getNumEdges();
		for (int edge = 0; edge < E; edge++) {
			int T = (*g)[source].getEdgeTarget(edge);
			if (T == target) {
				rates.push_back((*g)[source].getEdgeRate(edge));
				break;
			}
		}
	}
}


void getEndDistribution(Graph* g, int source, double* probs, int possible, std::map<int,int> toIndex, int* indices) {
	//get the probability of each end state for given graph g

	for (int i = 0; i < g->getN(); i++) {
		int E = (*g)[i].getNumEdges();
		int B = (*g)[i].getNumBackEdges();
		if (E == 0 && B != 0) {
			int inArray = toIndex.find(i)->second;
			probs[inArray] = (*g)[i].getProb();
		}
	}

	//if the node is an end node, this wont work. check seperately
	for (int i = 0; i < possible; i++) {
		if (indices[i] == source) {
			probs[i] = 1;
		}
	}

}

void findConditionalEnd(Graph* g, int source) {
	//for every node in g, find probability to end in each possible state
	//probabilities will be ordered in a vector, smallest to largest index

	//make arrays for storing info
	int possible = g->getNumEndStates();
	int* indices = new int[possible];
	double* probs = new double[possible];
	std::map<int, int> toIndex;
	std::vector<double> endDistr; 

	//do the worm state, standard graph. fill indices with end states.
	int count = 0;
	for (int i = 0; i < g->getN(); i++) {
		int E = (*g)[i].getNumEdges();
		int B = (*g)[i].getNumBackEdges();
		if (E == 0 && B != 0) {
			indices[count] = i; probs[count] = (*g)[i].getProb();
			toIndex.insert(std::make_pair(i, count));
			count++;
		}
	}
	for (int i = 0; i < possible; i++) {
		endDistr.push_back(probs[i]);
		probs[i] = 0;
	}
	(*g)[source].endDistr = endDistr;
	endDistr.clear();

	//make every subgraph, compute the new end distribution, fill it in. 
	for (int i = 0; i < g->getN(); i++) {
		if (i != source) {
			Graph* s = sourceSubgraph(g, i); 
			getEndDistribution(s, i, probs, possible, toIndex, indices);
			for (int j = 0; j < possible; j++) {
				endDistr.push_back(probs[j]);
				probs[j] = 0;
			}
			(*g)[i].endDistr = endDistr;
			endDistr.clear();
			delete s;
		}
	}

	//free memory
	delete []indices; delete []probs;
}











}