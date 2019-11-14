#include <math.h>
#include <stdio.h>
#include "nauty.h"
#include "database.h"
#include "adjacency.h"
#include "../defines.h"
#include <vector>
#include <eigen3/Eigen/Dense>
namespace bd{

class Database;

int toIndex(int r, int c, long m) {
  //map row and column number into index in 1d array. column indexed
  return m*c+r;
}

void index2ij(int index, int N, int& i, int& j) {
  // map a row index 1-d array into (i,j) 
  i = index % N;
  j = index / N;
}

void extractAM(int N, int state, int* AM, Database* db) {
	//extracts the adjacency matrix of state from the database

	//zero out the matrix from previous searches
	for (int i = 0; i < N*N; i++) AM[i] = 0;

	//fill the matrix with db data
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			AM[i*N+j] = (*db)[state].isInteracting(i,j,N);
		}
	}
}

void printAM(int N, int* AM) {
	//print out the adj matrix 
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf("%d ", AM[toIndex(N,i,j)]);
		}
		printf("\n");
	}
}

bool checkConnected(int* M, int N) {
	//check if subdiagonal sums to N-1
	int S = 0;
	for (int i = 0; i < N-1; i++) {
		S += M[(N+1)*i+1];
	}
	if (S == N-1) {
		return true;
	}
	else {
		return false;
	}
}

bool checkSame(int* M1, int* M2, int N) {
	//check if adjacency matrices M1 and M2 are the same
	for (int i = 0; i < N*N; i++) {
		if (M1[i]-M2[i] != 0) {
			return false;
		}
	}
	return true;
}

void getAdj(double* X, int N, int* M) {
	//get the adjacnecy matrix from a cluster X

	int i, j; 
	double* particles = new double[DIMENSION*N]; double* Z = new double[DIMENSION];
	c2p(X, particles, N);
	for (int index = 0; index < N*N; index++) {
		index2ij(index, N, i, j);
		if (j > i) {
			double d = euDist(particles, i, j, N, Z);
			if (d < BOND_CUTOFF) {
				M[index] = 1; M[toIndex(j,i,N)] = 1;
			}
		}
	}
	delete []particles; delete []Z;
}

void buildNautyGraph(int N, int M, int state, Database* db, graph* g) {
	//build nauty graph of given state

	//zero out any pre-existing graph
	EMPTYGRAPH(g, M, N);

	//add edges for every non-zero entry in adjacnecy matrix
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			if ((*db)[state].isInteracting(i,j,N)) {
				ADDONEEDGE(g, i, j, M);
			}
		}
	} 
}

bool checkIsomorphic(int N, int M, graph* g1, graph* g2) {
	//check if two canonically labeled graphs are isomorphic

	size_t k;
	for (k = 0; k < M*(size_t)N; ++k) {
			if (g1[k] != g2[k]) { //label is different, return 0
				return false;
			}
		}
		if (k == M*(size_t)N) { //graphs are isomorphic, return 1
			return true;
		}

}


void findIsomorphic(int N, int num_states, int state, Database* db, std::vector<int>& iso) {
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
	buildNautyGraph(N, M, state, db, g1);
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

void makeNM(int N, int* M, Eigen::VectorXd x, Eigen::MatrixXd& J, Eigen::VectorXd& F) {
	//make the matrix and vector to solve for Newtons method

	double XD, YD;
	int count = 0;

	//loop over and construct system
	for (int i = 0; i <N; i ++) {
		for (int j = i+1; j < N; j++) {
			if (M[toIndex(i,j,N)] == 1) {
				XD = x(2*i) - x(2*j);
				YD = x(2*i+1) - x(2*j+1);
				F(count) = XD*XD + YD*YD -1.0;
				J(count, 2*i) = 2*XD; J(count, 2*j) = -2*XD;
				J(count, 2*i+1) = 2*YD; J(count, 2*j+1) = -2*YD;
				count +=1;
			}
		}
	}
}

void refine(int N, double* X, int* M) {
	//refine a potentially unphysical state with NM - fill an adjacency matrix

	//set parameters to newtons method
	int max_iter = 50;
	double tol = 1e-8;

	int b = 0; //num bonds
	for (int i = 0; i < N*N; i++) b+=M[i];
	b /= 2;

	//initialize the matrix and vectors
	Eigen::MatrixXd J(b,2*N); 
	Eigen::VectorXd F(b); 
	Eigen::VectorXd dx(2*N);  
	Eigen::VectorXd x(2*N); 

	//initialize x. fill others with zeros.
	for (int i = 0; i < 2*N; i++) {
		x(i) = X[i];
	}
	dx.fill(0.0); F.fill(0.0); J.fill(0.0);

	//fill F and J. do solve with svd decomp.
	int iter;
	for (iter = 0; iter < max_iter; iter++) {
		makeNM(N, M, x, J, F);
		dx = J.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-F);
		x = x+dx;
		if (dx.norm() < tol) {
			break;
		}
	}

	//build adjacency matrix
	for (int i = 0; i < N*N; i++) M[i] = 0;
	double* XM = new double[2*N]; for (int i = 0; i < 2*N; i++) XM[i]=x(i);
	getAdj(XM, N, M);

	//free memory
	delete []XM;
}


}