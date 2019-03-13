#include <math.h>
#include "bDynamics.h"
namespace bd{


int toIndex(int r, int c, long m) {
  //map row and column number into index in 1d array. column indexed
  return m*c+r;
}

void index2ij(int index, int N, int& i, int& j) {
  // map a row index 1-d array into (i,j) 
  j = index % N;
  i = index / N;
}

int checkConnected(int* M, int N) {
	//check if subdiagonal sums to N-1
	int S = 0;
	for (int i = 0; i < N-1; i++) {
		S += M[(N+1)*i+1];
	}
	if (S == N-1) {
		return 1;
	}
	else {
		return 0;
	}
}

int checkSame(int* M1, int* M2, int N) {
	//check if adjacency matrices M1 and M2 are the same
	for (int i = 0; i < N*N; i++) {
		if (M1[i]-M2[i] != 0) {
			return 0;
		}
	}
	return 1;
}

void getAdj(double* X, int N, int*& M) {
	//get the adjacnecy matrix from a cluster X
	int i, j; double tol = 1e-5;
	double* particles = new double[2*N]; double* Z = new double[2];
	c2p(X, particles, N);
	for (int index = 0; index < N*N; index++) {
		index2ij(index, N, i, j);
		if (j > i) {
			double d = euDist(particles, i, j, N, Z);
			if (d < 1.1 + tol) {
				M[index] = d; M[toIndex(j,i,N)] = d;
			}
		}
	}
	delete []particles; delete []Z;
}

}