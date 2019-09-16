/* This code is to test the monte carlo on manifold sampler on a trimer system */ 

#include <cstdlib>
#include <stdio.h>
#include "bDynamics.h"
#include "sampling.h"
#include "../defines.h"

int main(int argc, char* argv[]) {

	//set parameters
	int N = 3;
	int num_samples = 1000;

	//set the initial and final state storage
	double* X0 = new double[DIMENSION*N];
	bd::setupChain(X0, N);

	//output the final state
	bd::printCluster(X0, N);

	//free memory
	delete []X0;

	//exit
	return 0;

}