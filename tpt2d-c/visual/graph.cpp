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


namespace bd{





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

void makeEdge(std::ofstream& out_str, int source, int target, double edgeWidth, double rate) {
	//draw an edge from source to target

	std::string s = std::to_string(rate);
	s.erase(s.end()-3,s.end());

	out_str << "\"" + std::to_string(source) + "\" -- \"" + 
	std::to_string(target) + "\" [penwidth = " + std::to_string(edgeWidth) + 
	", label = " + s + "]\n";
}

void graphP(std::ofstream& out_str, Database* db, int state, std::vector<Pair> P, int draw, int reduce, double normalizer, std::vector<int>& drawn) {
	//use the data in P for state to make nodes and edges

	//declare variables
	std::vector<int> rankVec; rankVec.clear();
	int index; double value;
	double mfpt; double S;
	std::vector<double> rates; rates.clear();
	std::vector<double> edgeWidth; edgeWidth.clear();

	//compute the rates for each state in the row
	mfpt = (*db)[state].getMFPT();
	S = (*db)[state].sumP();
	for (int i = 0; i < P.size(); i++) {
		index = P[i].index; value = P[i].value/S; //index and probability to go to i. 
		rates.push_back( (1/mfpt) * value);
		if (draw == 0) {
			edgeWidth.push_back(rates[i]/normalizer*30);
		}
		else if (draw == 1 && reduce == 0) {
			edgeWidth.push_back(rates[i]/normalizer*50);
		}
		else if (draw == 1 && reduce == 1) {
			edgeWidth.push_back(P[i].value/S*25);
		}
	}

	double tol = 0.1; //tolerance for printing an edge. 
	//print the node and edge to file
	for (int i = 0; i < P.size(); i++) {
		if (reduce == 0 || P[i].value/S > tol) {
			index = P[i].index;
			printCluster(out_str, index, draw);
			makeEdge(out_str, state, index, edgeWidth[i], rates[i]);
			drawn.push_back(index);
			rankVec.push_back(index); 
		}
	}
	sameRank(out_str, rankVec);
}

bool downGraph(std::ofstream& out_str, Database* db, int bonds, int draw, int reduce, std::vector<int>& drawn) {
	//search for states with given bonds and construct graph

	//search for states first and get the normalizing constant for rates
	bool flag = false; double normalizer = 0;
	for (int i = 0; i < db->getNumStates(); i++) {
		int B = (*db)[i].getBonds();
		if (B == bonds) { //state with B bonds found
			normalizer += 1 / (*db)[i].getMFPT();
			flag = true;
		}
	}

	//if there was a state at the next level, go through and graph them
	if (flag) {
		for (int i = 0; i < db->getNumStates(); i++) {
			int B = (*db)[i].getBonds();
			if (B == bonds) {
				//check if state has been drawn before edge making
				if (reduce == 0 || std::find(drawn.begin(), drawn.end(), i) !=drawn.end()) {
 					std::vector<Pair> P = (*db)[i].getP();
					graphP(out_str, db, i, P, draw, reduce, normalizer, drawn);
				}
			}
		}
	}
	return flag; //if flag is false, there are no states with B bonds, stop
}







void makeGraphViz(Database* db, int draw, int reduce) {
	//create a graphviz graph using mfpt data

	//get problem data
	int N = db->getN(); int ns = db->getNumStates();

	//declare variables
	int index; double value; 
	std::vector<int> rankVec;
	std::vector<int> drawn; 

	//declare a file to write output to
	std::string out;
	out = "N" + std::to_string(N) + "graphviz.txt";
	std::ofstream out_str(out);

	//write the graphviz header
	out_str << "graph mfptN" + std::to_string(N) + " {\n nodesep = 1.5; ranksep = 4; \n";
	out_str << "edge [ fontcolor=red, fontsize=36];\n";

	//set the minimum number of bonds - N-1
	int min_bond = N-1; int first = 0;

	//search the database for the worm state and create the node
	for (int i = 0; i < ns; i++) {
		if ((*db)[i].getBonds() == min_bond) {
			first = i; 
			break;
		}
	}
	printCluster(out_str, first, draw);

	//make a same rank statement for vertical alignment
	rankVec.push_back(first);
	sameRank(out_str, rankVec);
	rankVec.clear(); drawn.clear();

	//get the transition data and make new nodes and edges
	std::vector<Pair> P = (*db)[first].getP();
	graphP(out_str, db, first, P, draw, reduce, 1/(*db)[first].getMFPT(), drawn);

	//search the database for N+ bonded states and repeat
	for (int bond_num = N; bond_num < N*N; bond_num++) {
		bool stop = downGraph(out_str, db, bond_num, draw, reduce, drawn);
		if (!stop) {
			break;
		}
	}

	//end reached - put a curly to end
	out_str << "}";	
	

}

}