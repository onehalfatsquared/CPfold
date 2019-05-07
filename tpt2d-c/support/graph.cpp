#include "database.h"
#include "graph.h"
#include "pair.h"
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
			int source = (*g)[node].getEdgeSource(edge);
			if (inSub[source] == 0) {
				inSub[source] = 1;
				toVisit.push_back(source);
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
				v.prob = v.prob + (*g)[source].prob * prob;
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
					v.prob = v.prob + (*g)[node].prob * prob;
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

void printCluster(std::ofstream& out_str, int index, int draw) {
	//print a node creation in graphviz

	if (draw == 0) {
		out_str << "\"" + std::to_string(index) + "\" [label= \"" + std::to_string(index)
		+ "\" , shape=circle, width = 2, regular = 1, style = filled, fillcolor=white]; \n";
	}
	else if (draw == 1) {
		out_str << "\"" + std::to_string(index) + "\" [label=\"\", shape=circle, width = 1, regular = 1," +
		" style = filled, fillcolor=white," + "image=\"c" +std::to_string(index) + ".png\"]; \n";
	}
} 

void sameRank(std::ofstream& out_str, std::vector<int> states) {
	//print a rank = same statement for all elements of states

	out_str << "{rank = same; ";
	for (int i = 0; i < states.size(); i++) {
		out_str << "\"" + std::to_string(states[i]) + "\";";
	}
	out_str << "}\n";
}

void makeEdgeClean(std::ofstream& out_str, int source, int target, double edgeWidth) {
	//draw an edge from source to target - no labels

	out_str << "\"" + std::to_string(source) + "\" -- \"" + 
	std::to_string(target) + "\" [penwidth = " + std::to_string(edgeWidth) + "]\n";
}

void makeEdge(std::ofstream& out_str, int source, int target, double edgeWidth, double rate) {
	//draw an edge from source to target

	std::string s = std::to_string(rate);
	s.erase(s.end()-3,s.end());

	out_str << "\"" + std::to_string(source) + "\" -- \"" + 
	std::to_string(target) + "\" [penwidth = " + std::to_string(edgeWidth) + 
	", label = " + s + "]\n";
}

void printGraph(Graph* g, int source, int draw, int clean) {
	//construct a graphviz file to print the graph

	//keep track of drawn nodes and ranks
	std::vector<int> rankVec; std::deque<int> toVisit;
	int ns = g->getN();
	bool* drawn = new bool[ns]; for (int i = 0; i < ns; i++) drawn[i] = 0; 

	//declare a file to write output to
	std::string out;
	out = "graphviz.txt";
	std::ofstream out_str(out);

	//write the graphviz header
	out_str << "graph bd {\n nodesep = 1.5; ranksep = 4; \n";
	out_str << "edge [ fontcolor=red, fontsize=48];\n";

	//print the source node
	printCluster(out_str, source, draw);
	rankVec.push_back(source);
	sameRank(out_str, rankVec);
	rankVec.clear(); drawn[source] = 1;

	//loop over edges of source
	int E = (*g)[source].getNumEdges();
	double node_prob = (*g)[source].getProb();
	for (int edge = 0; edge < E; edge++) {
		int T = (*g)[source].getEdgeTarget(edge);
		double rate = (*g)[source].getEdgeRate(edge);
		double prob = (*g)[source].getEdgeProb(edge);
		if (!drawn[T]) {
			printCluster(out_str, T, draw);
			drawn[T] = 1; rankVec.push_back(T);
		}
		double edgeWidth = 40*node_prob*prob;
		if (clean == 0) {
			makeEdge(out_str, source, T, edgeWidth, rate);
		}
		else if (clean == 1) {
			makeEdgeClean(out_str, source, T, edgeWidth);
		}
		toVisit.push_back(T);
	}
	sameRank(out_str, rankVec); int rankNum = rankVec.size();
	rankVec.clear();

	//loop over the rest of the nodes
	int count = 0;
	while (!toVisit.empty()) {
		int node = toVisit[0];
		count++;
		int E = (*g)[node].getNumEdges();
		double node_prob = (*g)[node].getProb();
		for (int edge = 0; edge < E; edge++) {
			int T = (*g)[node].getEdgeTarget(edge);
			double rate = (*g)[node].getEdgeRate(edge);
			double prob = (*g)[node].getEdgeProb(edge);
			if (rate > 0.01) {
				if (!drawn[T]) {
					printCluster(out_str, T, draw);
					drawn[T] = 1; rankVec.push_back(T);
					toVisit.push_back(T);
				}
				double edgeWidth = 40*node_prob*prob;
				if (clean == 0) {
					makeEdge(out_str, node, T, edgeWidth, rate);
				}
				else if (clean == 1) {
					makeEdgeClean(out_str, node, T, edgeWidth);
				}
			}
		}
		toVisit.pop_front();
		if (count == rankNum) {
			count = 0; rankNum = rankVec.size();
			sameRank(out_str, rankVec); rankVec.clear();
		}
	}













	//end reached - put a curly to end
	out_str << "}";	
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

void printPath(std::vector<int> path, std::vector<double> val, std::string s) {
	//makes a graphviz file that prints the given path - with images, thats the point

	//declare a file to write output to
	std::string out;
	out = s + "graphviz.txt";
	std::ofstream out_str(out);

	//write the graphviz header
	out_str << "graph " + s + " {\n nodesep = 1.5; ranksep = 4; \n";
	out_str << "edge [ fontcolor=red, fontsize=48];\n";

	//draw all the nodes
	for (int i = 0; i < path.size(); i++) {
		int index = path[i];
		out_str << "\"" + std::to_string(index) + "\" [label=\"\", shape=circle, width = 1, regular = 1," +
			" style = filled, fillcolor=white," + "image=\"c" +std::to_string(index) + ".png\"]; \n";
	}

	//draw all the edges
	for (int i = 0; i < path.size()-1; i++) {
		int source = path[i]; int target = path[i+1];
		std::string s = std::to_string(val[i]);
		s.erase(s.end()-3,s.end());
		int edgeWidth = 5;
		out_str << "\"" + std::to_string(source) + "\" -- \"" + 
			std::to_string(target) + "\" [penwidth = " + std::to_string(edgeWidth) + 
			", label = " + s + "]\n";
	}

	//make all nodes have the same rank
	out_str << "{rank = same; ";
	for (int i = 0; i < path.size(); i++) {
		out_str << "\"" + std::to_string(path[i]) + "\";";
	}
	out_str << "}\n";

	//print ending curly
	out_str << "}";
}

void MPP(Graph* g, int source) {
	//determine the most probable path down the graph from given source

	//find the path by traveling from top to bottom
	std::vector<int> path; std::vector<double> prob;
	path.push_back(source);
	int E = (*g)[source].getNumEdges();
	while (E > 0) {
		findMaxProb(g, source, E, path, prob);
		E = (*g)[source].getNumEdges();
	}

	//print to user
	printf("Most Probable Path: ");
	for (int i = 0; i < path.size()-1; i++) {
		printf("%d, ", path[i]);
	}
	printf("%d\n", path[path.size()-1]);

	//make a graphviz file that shows the path
	printPath(path, prob, "MPP");

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

void QP(Graph* g, int source) {
	//determine the quickest folding path (largest rate)

	//get number of nodes
	int ns = g->getN();

	//make arrays that will store max sums and paths
	double* M = new double[ns];
	std::vector<std::vector<int>> paths(ns);
	for (int i = 0; i < ns; i++) M[i] = 0;
	M[source] = 0; paths[source].push_back(source);

	//compute a topological ordering of the graph
	int* T = new int[ns];
	tOrder(g, T, source);

	//loop over nodes, go down edges, check for largest path
	for (int i = 0; i < ns; i++) {
		int node = T[i]; 
		int E = (*g)[node].getNumEdges();
		for (int edge = 0; edge < E; edge++) {
			int target = (*g)[node].getEdgeTarget(edge);
			double rate = (*g)[node].getEdgeRate(edge);
			double pathTotal = M[node] + rate;
			if (pathTotal > M[target]) {
				M[target] = pathTotal;
				std::vector<int> prevPath = paths[node];
				prevPath.push_back(target);
				paths[target] = prevPath;
			}
		}
	}
	
	//find the max of the array M, get corresponding path
	int ind; double maxVal = 0;
	for (int i = 0; i < ns; i++) {
		if (M[i] > maxVal) {
			maxVal = M[i];
			ind = i;
		}
	}
	std::vector<int> maxPath = paths[ind];

	//print to user
	printf("Quickest Path: ");
	for (int i = 0; i < maxPath.size()-1; i++) {
		printf("%d, ", maxPath[i]);
	}
	printf("%d. ", maxPath[maxPath.size()-1]);
	printf("Total Rate = %f\n", M[ind]);

	//fill the rates on the path
	std::vector<double> rates;
	fillRates(g, maxPath, rates);

	//get graphviz code to show path
	printPath(maxPath, rates, "QP");

	//free memory
	delete []T; delete []M;
}

void QP(Graph* g, int source, int target) {
	//determine the quickest folding path (largest rate) from source to target

	//get number of nodes
	int ns = g->getN();

	//make arrays that will store max sums and paths
	double* M = new double[ns];
	std::vector<std::vector<int>> paths(ns);
	for (int i = 0; i < ns; i++) M[i] = 0;
	M[source] = 0; paths[source].push_back(source);

	//compute a topological ordering of the graph
	int* T = new int[ns];
	tOrder(g, T, source);

	//loop over nodes, go down edges, check for largest path
	for (int i = 0; i < ns; i++) {
		int node = T[i]; 
		int E = (*g)[node].getNumEdges();
		for (int edge = 0; edge < E; edge++) {
			int target = (*g)[node].getEdgeTarget(edge);
			double rate = (*g)[node].getEdgeRate(edge);
			double pathTotal = M[node] + rate;
			if (pathTotal > M[target]) {
				M[target] = pathTotal;
				std::vector<int> prevPath = paths[node];
				prevPath.push_back(target);
				paths[target] = prevPath;
			}
		}
	}
	
	//Get the desired target max
	int ind = target; double maxVal = M[target];
	std::vector<int> maxPath = paths[ind];

	//print to user
	printf("Quickest Path ending at %d: ", target);
	for (int i = 0; i < maxPath.size()-1; i++) {
		printf("%d, ", maxPath[i]);
	}
	printf("%d. ", maxPath[maxPath.size()-1]);
	printf("Total Rate = %f\n", M[ind]);

	//fill the rates on the path
	std::vector<double> rates;
	fillRates(g, maxPath, rates);

	//get graphviz code to show path
	printPath(maxPath, rates, "QP");

	//free memory
	delete []T; delete []M;
}










}