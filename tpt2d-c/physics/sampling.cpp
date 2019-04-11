#include <math.h>
#include <stdio.h>
#include <string.h>
#include "bDynamics.h"
#include "database.h"
#include <omp.h>
namespace bd{


double sampleSTD(double* X, int n) {
	//evaluate the sample standard deviation of data in X
	double mean = 0; double std = 0;
	for (int i = 0; i < n; i++) mean += X[i];
	mean /= n;
	for (int i = 0; i < n; i++) std += (X[i]-mean)*(X[i]-mean);
	std /= (n-1); 
	return sqrt(std);
}

void extractAM(int N, int state, int* AM, Database* db) {
	//extracts the adjacency matrix of state from the database
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			AM[i*N+j] = (*db)[state].isInteracting(i,j,N);
		}
	}
}

bool findMatrix(int* M, int* old, int old_bonds, int N, Database* db, int& timer, int& reset, int& reflect, 
																					int& new_state) {
	//find a state by transition matrix in the database

	for (int i = 0; i < (*db).getNumStates(); i++) {
		extractAM(N, i, old, db);
		int S = checkSame(old, M, N); 
		if (S == 1) {//found matrix
			int new_bonds = (*db)[i].getBonds();
			if (new_bonds == old_bonds + 1 || (old_bonds == 10 && new_bonds ==12 )) {//keep these
				new_state = i; reflect = 1;
			}
			else {// 2 states at once transition. just delete this sample
				reset = 1; timer = 0;
			}
			return 1; 
		}
	}
	return 0;
}


void checkState(double* X, int N, int state, int& new_state, Database* db, int& timer,
							 int& reset, int& reflect) {
	//check if state is different than prev step. if so, do various things

	//update timer
	timer +=1;

	//get adjacnecy matrix of the current state
	int* M = new int[N*N]; for (int i = 0; i < N*N; i++) M[i]=0;
	getAdj(X, N, M);

	//check if state is connected
	int C = checkConnected(M, N);
	if (C == 0) {//not connected
		reset = 1; timer -= 1;
		delete []M;
		return;
	}

	//compare old and new state, check if same
	int* old = new int[N*N]; for (int i = 0; i < N*N; i++) old[i]=0;
	int old_bonds = (*db)[state].getBonds();
	extractAM(N, state, old, db);
	int S = checkSame(old, M, N);

	if (S == 1) {//same matrix
		delete []M; delete []old;
		return;
	}
	else {//not the same, find matrix in database
		S = findMatrix(M, old, old_bonds, N, db, timer, reset, reflect, new_state);
		if (S == 0) {
			//most of these are unphysical 11 bond state - map to 12 bond with newton?
			for (int j = 0; j < 2*N; j++) printf("%f\n", X[j]);
			printf("State not found in database\n\n\n\n\n\n\n");
			reset = 1; timer = 0;
		}
		delete []M; delete []old;
		return;
	}
}


void runTrajectoryMFPT(double* X, Database* db, int state, int samples, int N, 
	double DT, int rho, double* E, double beta, int* P, int method, int& Num, 
																											int& Den, int* PM ) {
	//run the trajectory, update mfpt estimates

	//intiailize temp storage and set parameters
	double* temp = new double[2*N]; memcpy(temp, X, 2*N*sizeof(double));
	int reset; int reflect; int new_state = state; int hit = 0; int max_it = 10*samples;
	int timer = 0; 

	//solve sde and update
	for (int i = 0; i < max_it; i++) {
		reset = 0; reflect = 0;
		//solve SDE
		solveSDE(X, N, DT, rho, beta, E, P, method);
		//check if state changed
		checkState(X, N, state, new_state, db, timer, reset, reflect);
		if (reflect == 0 && reset == 0) {//no hit, proceed
			memcpy(temp, X, 2*N*sizeof(double)); // copy x to temp
		}
		else if (reflect == 1) {//hit new state, update estimates
			Den += timer; Num += timer*(timer+1)/2.0; 
			PM[new_state] += 1;
			memcpy(X, temp, 2*N*sizeof(double));//copy temp to x -> reset step
			timer = 0; hit +=1;
		}
		else {//chain broke, reset previous config
			memcpy(X, temp, 2*N*sizeof(double));//copy temp to x -> reset step
		}
		//if we reach desired number of samples, break
		if (hit == samples) {
			break;
		}
	}
	//free the temp memory
	delete []temp;
}


void equilibrate(double* X, Database* db, int state, int eq, int N, double DT, 
													int rho, double* E, double beta, int* P, int method) {
	//perform eq steps to equilibrate the trajectory. do not record data

	//initialize temp storage and set parameters;
	double* temp = new double[2*N]; memcpy(temp, X, 2*N*sizeof(double));
	int reset; int reflect; int new_state = state; int timer = 0;

	//solve the SDE and update
	for (int i = 0; i < eq; i++) {
		reset = 0; reflect = 0;
		//solve the sde
		solveSDE(X, N, DT, rho, beta, E, P, method);
		//check if state changed
		checkState(X, N, state, new_state, db, timer, reset, reflect);
		//if state changed, reflect back. otherwise continue
		if (reflect == 0 && reset == 0) {
			memcpy(temp, X, 2*N*sizeof(double)); // copy x to temp
		}
		else {
			memcpy(X, temp, 2*N*sizeof(double));//copy temp to x -> reset step
		}
	}
	//free the temp memory
	delete []temp;
}


void setupSimMFPT(int N, double Eh, int*& P, double*& E) {
	//initialize the interaction matrices - all interacting, strong bonds
	for (int i = 0; i < N*N; i++) {
		P[i] = 1; E[i] = Eh;
	}
}


void estimateMFPT(int N, int state, Database* db) {
	/*estimate mean first passage time starting in state and going to state with
	one additional bond. Uses parallel implementations of a single walker with
	long trajectory.*/

	//set parameters
	int rho = 40; double beta = 1; double DT = 0.01; int Kh = 1850;
	int method = 1; //solve SDEs with EM
	int samples = 150; //number of hits per walker for estimator
	int eq = 200; //number of steps to equilibrate for

	//quantities to update - old estimates
	int num_states = db->getNumStates();
	int num = (*db)[state].getNumerator(); 
	int den = (*db)[state].getDenominator(); 
	int* pm = new int[num_states];
	for (int i = 0; i < num_states; i++) {
		pm[i] = (*db)[state].getP(i);
	}

	//quantities to update - new estimates
	int NUM = 0; int DEN = 0;
	int* PM = new int[num_states];
	for (int i = 0; i < num_states; i++) {
		PM[i] = 0;
	}
	double mfpt = 0;
	double sigma = 0;

	//output start message
	printf("Beginning MFPT Estimator for state %d out of %d.\n", state, num_states);

	//setup simulation
	double Eh = stickyNewton(8, rho, Kh); //get energy corresponding to kappa
	//initialize interaction matrices
	int* P = new int[N*N]; double* E = new double[N*N];
	setupSimMFPT(N, Eh, P, E);

	//store mfpt estimates on each thread to get standard deviation
	double* mfptSamples; int num_threads;

	//open parallel region
	#pragma omp parallel reduction(+:NUM, DEN, PM[:num_states])
	{

	//initialize samples
	num_threads = omp_get_num_threads();
	mfptSamples = new double[num_threads];

	//get starting structures
	Cluster c = (*db)[state].getRandomIC();

	//cluster structs to arrays
	double* X = new double[2*N];
	c.makeArray2d(X, N);

	//equilibrate the trajectories
	equilibrate(X, db, state, eq, N, DT, rho, E, beta, P, method);

	//run BD
	runTrajectoryMFPT(X, db, state, samples, N, DT, rho, E, beta, P, method, NUM, DEN, PM );

	//store samples
	if (DEN != 0) {
		mfptSamples[omp_get_thread_num()] = (NUM * DT) / DEN;
	}
	else {
		mfptSamples[omp_get_thread_num()] = 0;
	}

	//free cluster memory
	delete []X;

	//end parallel region
	}

	//combine estimates - if any samples were found (11->12 state)
	if (DEN != 0) {
		num += NUM; den += DEN;
		for (int i = 0; i < num_states; i++) {
			pm[i] += PM[i];
		}
		mfpt = (num * DT) / den;
		sigma = sampleSTD(mfptSamples, num_threads);

		//update database
		(*db)[state].mfpt = mfpt;
		(*db)[state].num = num;
		(*db)[state].denom = den;
		for (int i = 0; i < num_states; i++) {
			(*db)[state].P[i] = pm[i];
		}
		(*db)[state].sigma = sigma;
	}

	//print out final estimates - debug
	
	printf("%d, %d, \n", num, den);
	for (int i = 0; i < num_states; i++) {
		if (pm[i]>0)
			printf("%d, %d\n", pm[i], i);
	}
	printf("Estimate = %f +- %f\n", mfpt, sigma);
	for (int i = 0; i < num_threads; i++) printf("MFPT estimate %d = %f\n", i, mfptSamples[i]);
	


	//free memory
	delete []E; delete []P; delete []PM; delete []pm; delete []mfptSamples;
}


void estimatePartitionFn(int N, int state, Database* db) {
	/*estimate mean first passage time starting in state and going to state with
	one additional bond. Uses parallel implementations of a single walker with
	long trajectory.*/

	//set parameters
	int rho = 40; double beta = 1; double DT = 0.01; int Kh = 1850;
	int method = 1; //solve SDEs with EM
	int samples = 1000; //number of hits per walker for estimator
	int eq = 200; //number of steps to equilibrate for




}





}