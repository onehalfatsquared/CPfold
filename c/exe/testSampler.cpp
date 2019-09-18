/* This code is to test the monte carlo on manifold sampler on a trimer system */ 

#include <cstdlib>
#include <stdio.h>
#include "bDynamics.h"
#include "sampling.h"
#include "../defines.h"

int main(int argc, char* argv[]) {

	//set parameters
	int N = 3;
	int num_samples = 100000;
	int kappa = 1.5;

	//set the initial and final state storage
	double* X0 = new double[DIMENSION*N];
	bd::setupChain(X0, N);

	//test
	int* M = new int[N*N];
	for (int i = 0; i < N*N; i++) M[i] = 0;
	M[1] = 1; M[3] = 1; M[5] = 1; M[7] = 1;
	mcm::runTest(N,X0,M,num_samples);

	//output the final state
	bd::printCluster(X0, N);

	//free memory
	delete []X0; delete []M;

	//exit
	return 0;

}