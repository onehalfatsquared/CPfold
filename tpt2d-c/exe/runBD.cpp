#include <cstdlib>
#include <stdio.h>
#include "bDynamics.h"

int main(int argc, char* argv[]) {

	//handle input
	if (argc != 3) {
		fprintf(stderr, "Usage: %s <N> <T>", argv[0]);
		return 1;
	}
	int N = atoi(argv[1]);
	double T = atof(argv[2]);

	//set parameters
	int rho = 40;
	double beta = 1;
	int method = 1;

	//initialize interaction matrices
	int* P = new int[N*N];
	double* E = new double[N*N];
	for (int i = 0; i < N*N; i++) {
		P[i] = 1; E[i] = 10;
	}

	//set the initial and final state storage
	double* X0 = new double[2*N];
	for (int i = 0; i < 2*N; i++) {
		if (i % 2 == 0) {
			X0[i] = i/2+1;
		}
		else{
			X0[i] = 0;
		}
	}


	bd::solveSDE(X0, N, T, rho, beta, E, P, method);

	for (int i = 0; i < 2*N; i++) printf("%f\n",X0[i]);


	delete []P;
	delete []E;
	delete []X0;
	return 1;

}
