#include <cstdlib>
#include <stdio.h>
#include "bDynamics.h"
#include "../defines.h"


int main(int argc, char* argv[]) {

	//handle input
	if (argc != 3) {
		fprintf(stderr, "Usage: %s <Num Particles> <Final Time>", argv[0]);
		return 1;
	}
	int N = atof(argv[1]);
	double T = atof(argv[2]);

	//set parameters
	int rho = RANGE;
	double beta = BETA;
	int method = 1;

	//initialize interaction matrices
	int* types = new int[N];
	int numTypes = bd::readDesignFile(N, types);
	int numInteractions = numTypes*(numTypes+1)/2;
	double* kappa = new double[numInteractions];
	bd::readKappaFile(numInteractions, kappa);
	std::map<std::pair<int,int>, double> kmap; 
	bd::makeKappaMap(numTypes, kappa, kmap);
	int* P = new int[N*N];
	double* E = new double[N*N];
	bd::fillP(N, types, P, E, kmap);

	//set the initial and final state storage
	double* X0 = new double[DIMENSION*N];
	bd::setupChain(X0, N);

	//set potential type
	int pot = 0; //0 morse, 1 lj

	//run bd
	bd::solveSDE(X0, N, T, rho, beta, E, P, method, pot);

	//output the final state
	bd::printCluster(X0, N);

	//free memory
	delete []P; delete []E; delete []X0; delete []kappa; delete []types;

	//exit
	return 0;

}
