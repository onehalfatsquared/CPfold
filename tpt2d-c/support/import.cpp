#include "database.h"
#include "nauty.h"
#include "adjacency.h"
#include <algorithm>
//#include <../Eigen/Dense>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
namespace bd {

void buildNautyGraph(int N, int M, int* AM, graph* g) {
	//build nauty graph of state with adj matrix M

	//zero out any pre-existing graph
	EMPTYGRAPH(g, M, N);

	//add edges for every non-zero entry in adjacnecy matrix
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			if (AM[j*N+i]) {
				ADDONEEDGE(g, i, j, M);
			}
		}
	} 
}

void findIsomorphic(int N, int num_states, int* AM, Database* db, std::vector<int>& iso) {
	//finds all states isomorphic to given state in database. stored to iso.

	//set keywords for nauty
	int M = SETWORDSNEEDED(N);

	//initialize the nauty graphs and parameters
	graph g1[N*M]; graph g2[N*M];   //the starting graphs
	graph cg1[N*M]; graph cg2[N*M]; //graphs with canonical labeling
	int lab1[N], lab2[N];           //label arrays
	int ptn[N], orbits[N];          //needed to call functions
	static DEFAULTOPTIONS_GRAPH(options); //defualt options
	statsblk stats;                 //nauty statistics of graph
	options.getcanon = TRUE;        //get canonical labeling
	options.defaultptn = TRUE;      //ignore any coloring of graph

	//build nauty graph of the state, get canonical labeling
	EMPTYGRAPH(cg1, M, N);
	buildNautyGraph(N, M, AM, g1);
	densenauty(g1, lab1, ptn, orbits, &options, &stats, M, N, cg1);

	//loop over all states, check if isomorphic
	for (int i = 0; i < num_states; i++) {
		//create graph
		EMPTYGRAPH(cg2, M, N);
		buildNautyGraph(N, M, i, db, g2);
		densenauty(g2, lab2, ptn, orbits, &options, &stats, M, N, cg2);

		/* Compare canonically labelled graphs */
		if (checkIsomorphic(N, M, cg1, cg2) == 1) { //graphs are isomorphic, add to iso
			iso.push_back(i);
		}
	} 
}

void structureMap(std::string& filename, Database* db, int* M2A, int& NC) {
	//read in data from an xyz file, create mapping from outside structures to mine

	std::ifstream in_str(filename);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		//return NULL;
	}

	//get the number of particles in the system
	int N;
	in_str >> N;
	int ns = db->getNumStates(); 

	//declare storage for cluster (x1,y1,x2,y2,...)
	double* X = new double[2*N]; 

	//declare storage for adjacency matrix
	int* M = new int[N*N];

	//declare a vector to store isoorphic states
	std::vector<int> iso;

	//declare storage for each line of file, and buffer to get values
	std::string line;
	char buffer[100];

	//counter to stop when cluster is filled
	int count = 0;
	int clust = 0;


	while (std::getline(in_str, line)) {
		//if line begins in "C", get line. extract coordinates
		if (line[0] == 'C') {
			//std::cout << line << "\n";
			strcpy(buffer,line.c_str());
			char* val = strtok(buffer, " \n");
			val = strtok(NULL, " \n");
			X[count] = atof(val);
			//printf("%f\n", X[count]);
			count++;
			val = strtok(NULL, " \n");
			X[count] = atof(val);
			//printf("%f\n", X[count]);
			count++;
		}

		if (count == 2*N) { //cluster is filled with coordinates
			//get the adjacency matrix for the cluster 
			//for (int i = 0; i < 2*N; i++) printf("aaaaa %f\n",X[i]);

			//do refining


			for (int i = 0; i < N*N; i++) M[i] = 0;
			getAdj(X, N, M);
			findIsomorphic(N, ns, M, db, iso);
			for (int i = 0; i < iso.size(); i++) {
				printf("Andreas Cluster %d is isomorphic to my cluster %d\n", clust, iso[i]);
				M2A[iso[i]] = clust;
			}
			if (iso.size() == 0) {
				printf("Andreas cluster %d has no match in my set.\n", clust);
			}
			iso.clear();
			clust++;
			count = 0;
		}
	}

	NC = clust;
	delete []X; delete []M;
}

void printTM(Database* db, int* M2A, int NC) {
	//print out the transition matrix to a file

	double tol = 0.001; //prob tolerance for considering a state to have transitions
	double* TM = new double[NC*NC];
	for (int i = 0; i < NC*NC; i++) {
		TM[i] = 0;
	}

	int ns = db->getNumStates();

	for (int i = 0; i < ns; i++) {
		std::vector<Pair> P = (*db)[i].getP();
		double S = (*db)[i].sumP();
		for (int entry = 0; entry < P.size(); entry++) {
			int index = P[entry].index; double val = P[entry].value;
			if (M2A[i] != -1 && M2A[index]!= -1 && val/S > tol) {
				TM[toIndex(M2A[i],M2A[index],NC)] = val/S; 
			}
		}
	}

	//print 
	std::string out = "one_step.txt";
	std::ofstream out_str(out);
	
	for (int i = 0; i < NC; i++) {
		for (int j = 0; j < NC; j++) {
			out_str << TM[toIndex(i,j,NC)] << ' ';
		}
		out_str << '\n';
	}

	delete [] TM;
}

void printTM(Database* db) {
	//print out the transition matrix to a file

	double tol = 0.001; //prob tolerance for considering a state to have transitions
	int ns = db->getNumStates();
	double* TM = new double[ns*ns];
	for (int i = 0; i < ns*ns; i++) {
		TM[i] = 0;
	}


	for (int i = 0; i < ns; i++) {
		std::vector<Pair> P = (*db)[i].getP();
		double S = (*db)[i].sumP();
		for (int entry = 0; entry < P.size(); entry++) {
			int index = P[entry].index; double val = P[entry].value;
			if (val/S > tol) {
				TM[toIndex(i,index,ns)] = val/S; 
			}
		}
	}

	//print 
	std::string out = "one_step_noMap.txt";
	std::ofstream out_str(out);
	
	for (int i = 0; i < ns; i++) {
		for (int j = 0; j < ns; j++) {
			out_str << TM[toIndex(i,j,ns)] << ' ';
		}
		out_str << '\n';
	}

	delete [] TM;
}












}