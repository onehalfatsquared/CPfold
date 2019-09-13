#include <cstdlib>
#include <stdio.h>
#include "bDynamics.h"
#include "../defines.h"

int main(int argc, char* argv[]) {

	//handle input
	if (argc != 2) {
		fprintf(stderr, "Usage: %s <Final Time>", argv[0]);
		return 1;
	}
	double T = atof(argv[1]);
	int N = PARTICLES;

	//set parameters
	int rho = RANGE;
	double beta = 1;
	int method = 1;

	//initialize interaction matrices
	int* P = new int[N*N];
	double* E = new double[N*N];
	for (int i = 0; i < N*N; i++) {
		P[i] = 1; E[i] = 10;
	}

	//set the initial and final state storage
	double* X0 = new double[DIMENSION*N];
	bd::setupChain(X0, N);

	//set potential type
	int pot = POTENTIAL; //0 morse, 1 lj

	//run bd
	bd::solveSDE(X0, N, T, rho, beta, E, P, method, pot);

	//output the final state
	for (int i = 0; i < DIMENSION*N; i++) printf("%f\n",X0[i]);

	//free memory
	delete []P; delete []E; delete []X0;

	//exit
	return 0;

}
