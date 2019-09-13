#include <random>
#include <math.h>
#include <chrono>
#include "bDynamics.h"
#include "../defines.h"
namespace bd{




void EM(double* X0, int N, int Nt, double k, 
					int rho, double* E, int* P, double beta, int pot) {
	//apply the EM method to solve the SDE
	//initialize particle and gradient storage
	double* g = new double[DIMENSION*N];
	double* particles = new double[DIMENSION*N];

	//initialize the random number generator - Normal(0,1) //fix seed
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<double> distribution(0.0,1.0);

	//apply the EM scheme
	for (int i = 0; i < Nt; i++) {
		c2p(X0, particles, N);
		if (pot == 0) {//use morse potential
			morseGrad(particles, rho, E, N, P, g);
		}
		else if (pot == 1) {//use lennard jones potential
			ljGrad(particles, rho, E, N, P, g);
		}
		for (int j = 0; j < DIMENSION*N; j++) {
			X0[j] += -g[j]*k + sqrt(2.0*k/beta)*distribution(generator);
		}
	}

	//free the memory
	delete []g; delete []particles;
}






void solveSDE(double* X0, int N, double T, int rho, double beta,
												 double* E, int* P, int method, int pot) {
	if (method == 1) {
		//set time step
		double k = EULER_TS; int Nt = T/k; 
		//solve the sde
		EM(X0,  N, Nt, k, rho, E, P, beta, pot);
	}
}

void setupChain(double* X, int N) {
	//construct a linear chain of particle positions
	for (int i = 0; i < DIMENSION*N; i++) {
		if (i % DIMENSION == 0) {
			X[i] = i/DIMENSION+1;
		}
		else{
			X[i] = 0;
		}
	}
}

}