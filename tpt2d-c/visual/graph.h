#pragma once
#include <vector>
#include <ios>
namespace bd { 
class Database; 


void makeGraphViz(Database* db, int draw, int reduce);
void sameRank(std::ofstream& out_str, std::vector<int> states);
void makeEdge(std::ofstream& out_str, int source, int target, double edgeWidth, double rate);
void printCluster(std::ofstream& out_str, int index, int draw);
void graphP(std::ofstream& out_str, Database* db, int state, 
							std::vector<Pair> P, int draw, int reduce, double normalizer, std::vector<int>&);
bool downGraph(std::ofstream& out_str, Database* db, int bonds, int draw,int reduce, std::vector<int>&);





}