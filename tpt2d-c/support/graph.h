#pragma once
#include <vector>
#include <ios>

/* Graph class to traverse and compute paths.
	Edge struct to store edge info:
		target - what directed edge points to
		rate - rate of going to target
		prob - probability to go to target

	Vertex class to store nodes in the graph
		index - label for node
		prob - probability to reach this state during folding
		edges - a vector of edge structs

*/

namespace bd { 

class Database; 

class Graph;

struct Edge {
	//null constructor
	Edge() {target = rate = prob = 0;}

	//define members of edge
	int target; double rate, prob;

	//constructor
	Edge(int target_, double rate_, double prob_) : target(target_), rate(rate_), prob(prob_) {}
};

class Vertex {
	public:
		Vertex();
		~Vertex();
		friend Graph* makeGraph(Database* db);
		friend void updateGraph(int node, Database* db, Graph* g, bool* present, int& count);

		//acessor functions
		int getIndex() const {return index;}
		double getProb() const {return prob;}
		int getNumEdges() const {return edges.size();}
		int getEdgeTarget(int i) const {return edges[i].target;}
		double getEdgeRate(int i) const {return edges[i].rate;}
		double getEdgeProb(int i) const {return edges[i].prob;}

	private:
		int index;
		double prob;
		std::vector<Edge> edges;

	//copy constructors - restricts compiling when user tries to copy a vertex
		Vertex(const Vertex&) {
			throw 1;
		}
		Vertex& operator=(const Vertex&) {
			throw 1;
		}


};

class Graph {
	public:
		Graph(int N_);
		~Graph();
		Vertex& operator[](int index) {return vertices[index];}
		const Vertex& operator[](int index) const {return vertices[index];}
		friend void updateGraph(int node, Database* db, Graph* g, bool* present, int& count);

		//accessor functions
		int getN() const {return N;}

	private:
		int N; Vertex* vertices;


	//copy constructors - restricts compiling when user tries to copy a graph
		Graph(const Graph&) {
			throw 1;
		}
		Graph& operator=(const Graph&) {
			throw 1;
		}
};


Graph* makeGraph(Database* db);
void endDistribution(Graph* g);
void MPP(Graph* g, int source);
int findMaxProb(Graph* g, int source, int E);
void QP(Graph* g, int source);








}