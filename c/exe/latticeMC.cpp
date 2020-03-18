#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "latticeP.h"
#include "../defines.h"

/* Run an mcmc simulation of a lattice protein  */

int main(int argc, char* argv[]) {

	//handle input
	if (argc != 3) {
		fprintf(stderr, "Usage: <Num particles> <useFile>  %s\n", argv[0]);
		return 1;
	}
	int N = atoi(argv[1]);
	bool useFile = atoi(argv[2]);

	lattice::runMCMC(N, useFile);

	return 0;
}