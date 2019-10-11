#pragma once
#include <vector>
#include <ios>
namespace bd { 
class Database; 
class Graph;


void printGraph(Graph* g, int source, int draw, int clean, int reduce) ;
void makeEdge(std::ofstream& out_str, int source, int target, double edgeWidth, double rate) ;
void makeEdgeClean(std::ofstream& out_str, int source, int target, double edgeWidth);
void sameRank(std::ofstream& out_str, std::vector<int> states);
void printCluster(std::ofstream& out_str, int index, int draw);

void printCluster(std::ofstream& out_str, int index, std::vector<double> end);
void printGraphEndDistribution(Graph* g, int source, int reduce);
//make a graphviz file with the given path
void printPath(std::vector<int> path, std::vector<double> val, std::string s);
//print the most probable path
void MPP(Graph* g, int source);
//print the quickest folding path
void QP(Graph* g, int source);
//print the quickest folding path ending at target
void QP(Graph* g, int source, int target);

//functions for reversible process - directed graphs
void printGraphRev(Graph* g, int source, double* F, double* flux, int draw, int clean, 
									 int reduce);
void makeEdgeRev(std::ofstream& out_str, int source, int target, double edgeWidth, 
								 double rate, std::string color);
void makeEdgeCleanRev(std::ofstream& out_str, int source, int target, double edgeWidth,
											std::string color);
void printClusterRev(std::ofstream& out_str, int index, double pw);

}
