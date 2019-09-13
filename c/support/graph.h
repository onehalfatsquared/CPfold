#pragma once
#include <vector>
#include <map>
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
		back_edges - a vector if ints corresponding to which nodes come into this node

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
		friend Graph* targetSubgraph(Graph* g, int source, int target);
		friend Graph* sourceSubgraph(Graph* g, int source);
		friend void findConditionalEnd(Graph* g, int source);

		//end state distributions
		std::vector<double> endDistr;

		//acessor functions
		int getIndex() const {return index;}
		double getProb() const {return prob;}
		int getNumEdges() const {return edges.size();}
		int getNumBackEdges() const {return back_edges.size();}
		int getEdgeTarget(int i) const {return edges[i].target;}
		int getEdgeSource(int i) const {return back_edges[i];}
		double getEdgeRate(int i) const {return edges[i].rate;}
		double getEdgeProb(int i) const {return edges[i].prob;}

	private:
		int index;
		double prob;
		std::vector<Edge> edges;
		std::vector<int> back_edges; 

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
		friend Graph* makeGraph(Database* db);

		//accessor functions
		int getN() const {return N;}
		int getNumEndStates() const {return end;}

	private:
		int N; int end;
		Vertex* vertices;


	//copy constructors - restricts compiling when user tries to copy a graph
		Graph(const Graph&) {
			throw 1;
		}
		Graph& operator=(const Graph&) {
			throw 1;
		}
};

//make the graph structure from the database
Graph* makeGraph(Database* db);
//make a subgraph of all states that can reach target
Graph* targetSubgraph(Graph* g, int source, int target);
//make a subgraph of all states after source
Graph* sourceSubgraph(Graph* g, int source);
//get the probability distribution of end states
void endDistribution(Graph* g);
//find the most probable path
void MPP(Graph* g, int source);
//find the outgoing edge with the largest value
void findMaxProb(Graph* g, int& source, int E, std::vector<int>& path, 
																					std::vector<double>& prob);
//make a graphviz file with the given path
void printPath(std::vector<int> path, std::vector<double> val, std::string s);
//find the quickest folding path
void QP(Graph* g, int source);
//find the quikcest folding path ending at target
void QP(Graph* g, int source, int target);
//construct a topological ordering of the graph
void tOrder(Graph* g, int* T, int source);
//given a path, fill in the rates for each step of the path
void fillRates(Graph* g, std::vector<int> path, std::vector<double>& rates);

void printGraph(Graph* g, int source, int draw, int clean, int reduce) ;
void makeEdge(std::ofstream& out_str, int source, int target, double edgeWidth, double rate) ;
void makeEdgeClean(std::ofstream& out_str, int source, int target, double edgeWidth);
void sameRank(std::ofstream& out_str, std::vector<int> states);
void printCluster(std::ofstream& out_str, int index, int draw);

void getEndDistribution(Graph* g, int source, double* probs, int possible, std::map<int,int> toIndex, int* indices);
void findConditionalEnd(Graph* g, int source);
void printCluster(std::ofstream& out_str, int index, std::vector<double> end);
void printGraphEndDistribution(Graph* g, int source, int reduce);








}