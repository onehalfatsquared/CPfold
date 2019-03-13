#include <math.h>
#include <stdio.h>
#include <string.h>
#include "bDynamics.h"
namespace bd{





void checkState(double* X, int N, int state, int new_state, ??? db, int timer,
							 int reset, int reflect) {
	//check if state is different than prev step. if so, do various things

	//update timer
	timer +=1;

	//get adjacnecy matrix of the current state
	int* M = new int[N*N]; for (int i = 0; i < N*N; i++) M[i]=0;
	getAdj(double* X, int N, int*& M);

	//check if state connected
	int C = checkConnected(M, N);
	if (C == 0) {//not connected
		reset = 1; timer -= 1;
		delete []M;
		return;
	}

	//compare old and new state, check if same
	int* old = ?; //todo
	int S = checkSame(old, M, N);
	if (S == 1) {//same matrix
		delete []M; // delete old as well? todo
		return;
	}
	else {//not the same, find matrix in database
		//todo

	}
}


void runTrajectory(double* X, int state, int samples, int N, double DT, int rho, double* E, 
					double beta, int* P, int method, double Num, double Den, double* PM ) {
	//run the trajectory, record data.
	double* temp = new double[2*N]; memcpy(temp, X, 2*N*sizeof(double));
	int reset; int reflect; int new_state = state; int hit = 0; int max_it = 1e6;
	for (int i = 0; i < max_it; i++) {
		reset = 0; reflect = 0;
		solveSDE(X, N, DT, rho, beta, E, P, method);
		//check if state changed
		checkState(X, N, state, new_state, ?, timer, reset, reflect);
		if (reflect == 0 && reset == 0) {//no hit, proceed
			memcpy(temp, X, 2*N*sizeof(double)); // copy x to temp
		}
		else if (reflect == 1) {//hit state, update estimates
			Den += timer; Num += timer*(timer+1)/2.0; timer = 0;
			PM[new_state] += 1;
			memcpy(X, temp, 2*N*sizeof(double));//copy temp to x -> reset step
		}
		else {//chain broke, reset previous config
			memcpy(X, temp, 2*N*sizeof(double));//copy temp to x -> reset step
		}
	}
	delete []temp;
}






void equilibrate(double* X, int state, int eq, int N, double DT, int rho, double* E, 
										double beta, int* P, int method) {
	//perform eq steps to equilibrate the trajectory. do not record data
	double* temp = new double[2*N]; memcpy(temp, X, 2*N*sizeof(double));
	int reset; int reflect; int new_state = state;
	for (int i = 0; i < eq; i++) {
		reset = 0; reflect = 0;
		solveSDE(X, N, DT, rho, beta, E, P, method);
		//check if state changed
		checkState(X, N, state, new_state, ?, timer, reset, reflect);
		if (reflect == 0 && reset == 0) {
			memcpy(temp, X, 2*N*sizeof(double)); // copy x to temp
		}
		else {
			memcpy(X, temp, 2*N*sizeof(double));//copy temp to x -> reset step
		}
	}
	delete []temp;
}



void setupSim(int N, double Eh, int*& P, double*& E) {
	//initialize the interaction matrices - all interacting, strong bonds
	for (int i = 0; i < N*N; i++) {
		P[i] = 1; E[i] = Eh;
	}
}


void estimateMFPT(int N, int state) {
	/*estimate mean first passage time starting in state and going to state with
	one additional bond. Uses parallel implementations of a single walker with
	long trajectory.*/

	//set parameters
	int rho = 40; double beta = 1; double DT = 0.01; int Kh = 1850;
	int method = 1; int timer = 0; 
	int samples = 2000; int eq = 200; 

	//quantities to update
	double Num; double Den; double* PM;

	//import structure
	//todo

	//setup simulation
	double Eh = stickyNewton(8, rho, Kh);
	//initialize interaction matrices
	int* P = new int[N*N]; double* E = new double[N*N];
	setupSim(N, Eh, P, E);

	//get starting structures
	//todo

	//equilibrate the trajectories
	//fn

	//run BD
	




	//free memory
	delete []E; delete []P;


}











}