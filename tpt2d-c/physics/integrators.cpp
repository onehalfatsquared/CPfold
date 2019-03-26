#include <random>
#include <math.h>
#include <chrono>
#include "bDynamics.h"
namespace bd{




void EM(double* X0, int N, int Nt, double k, 
					int rho, double* E, int* P, double beta) {
	//apply the EM method to solve the SDE
	//initialize particle and gradient storage
	double* g = new double[2*N];
	double* particles = new double[2*N];

	//initialize the random number generator - Normal(0,1) //fix seed
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<double> distribution(0.0,1.0);

	//apply the EM scheme
	for (int i = 0; i < Nt; i++) {
		c2p(X0, particles, N);
		morseGrad(particles, rho, E, N, P, g);
		for (int j = 0; j < 2*N; j++) {
			X0[j] += -g[j]*k + sqrt(2.0*k/beta)*distribution(generator);
		}
	}

	//free the memory
	delete []g; delete []particles;
}






void solveSDE(double* X0, int N, double T, int rho, double beta,
												 double* E, int* P, int method) {
	if (method == 1) {
		//set time step
		double k = 5e-6; int Nt = T/k; 
		//solve the sde
		EM(X0,  N, Nt, k, rho, E, P, beta);
	}
}

}