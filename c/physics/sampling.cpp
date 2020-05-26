#include <math.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <complex>
#include <cmath>
#include "bDynamics.h"
#include "sampling.h"
#include "database.h"
#include "../defines.h"
#include <omp.h>

namespace bd{


/******************************************************************/
/**************** Functions to build Database *********************/
/******************************************************************/

void buildEmptyDB(int N) {
	//construct a new database file that only has the linear state in it

	//make the db with only 1 state - linear chain
	Database* db = new Database(N, 1);

	//create reference to that state
	State& s = (*db)[0];

	//fill in the state info

	//adjacency matrix
	s.am = new bool[N*N];
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (j == i+1 || j == i-1) {
				s.am[toIndex(i,j,N)] = 1;
				s.am[toIndex(j,i,N)] = 1;
			}
			else {
				s.am[toIndex(i,j,N)] = 0;
				s.am[toIndex(j,i,N)] = 0;
			}
		}
	}

	//bonds
	s.bond = N-1;

	//coordinates
	s.num_coords = 1;
	s.coordinates = new Cluster[s.num_coords];
	for (int i = 0; i < s.num_coords; i++) {
		s.coordinates[i].setNumPoints(N);
		for (int j = 0; j < N; j++) {
#if (DIMENSION == 2)
			double x = j;
			double y = 0.0;
			s.coordinates[i][j] = Point(x,y);
#elif (DIMENSION == 3)
			double x = j;
			double y = 0.0;
			double z = 0.0;
			s.coordinates[i][j] = Point(x,y,z);
#endif
		}
	}

	s.num = s.denom = 0;
	s.num_neighbors = 0;
	s.freq = s.mfpt = s.sigma = 0.0;

	//print out the new database
	std::string out = "N" + std::to_string(N) + "DB.txt";
	std::ofstream out_str(out);
	out_str << *db;
	delete db;
}

bool checkSame(int N, int* AM, State& s) {
	//check if states are same by adjacency matrix

	for (int i = 0; i < N; i++) {
		for (int j = i+2; j < N; j++) {
			if (s.isInteracting(i,j,N) != AM[toIndex(i,j,N)]) {
				return false;
			}
		}
	}

	return true;
}

int searchDB(int N, Database* db, std::vector<State> new_states, int* AM) {
	//check if the current adj matrix is in db. if yes, return state #. if not, return -1.

	int num_states = db->getNumStates();

	for (int state = 0; state < num_states; state++) {
		State& s = (*db)[state];
		bool same = checkSame(N, AM, s);
		if (same)
			return state;
	}

	//std::cout << new_states.size() << "\n\n\n\n\n";
	for (int state = 0; state < new_states.size(); state++) {
		State& s = new_states[state];
		bool same = checkSame(N, AM, s);
		if (same)
			return state+num_states;
	}

	return -1;
}

void addState(int N, double* X, int* AM, std::vector<State>& new_states) {
	//add a new state to a vector

	State s = State();

	int b = 0;
	s.N = N;

	s.am = new bool[N*N];
	//construct AM
	for (int i = 0; i < N*N; i++) {
		s.am[i] = AM[i];
		if (AM[i] == 1) {
			b++;
			//std::cout <<"hello" << s.isInteracting(i / N, i % N) << "\n";
		}
	}

	//add num bonds
	b /= 2;
	s.bond = b;

	//add configuration
	int coord = 1;
	s.num_coords = coord;
	s.coordinates = new Cluster[coord];
	for (int i = 0; i < coord; i++) {
		s.coordinates[i].setNumPoints(N);
		for (int j = 0; j < N; j++) {
#if (DIMENSION == 2)
			double x = X[DIMENSION*j];
			double y = X[DIMENSION*j+1];
			s.coordinates[i][j] = Point(x,y);
#elif (DIMENSION == 3)
			double x = X[DIMENSION*j];
			double y = X[DIMENSION*j+1];
			double z = X[DIMENSION*j+2];
			s.coordinates[i][j] = Point(x,y,z);
#endif
		}
	}

	//push to vector
	new_states.push_back(s);
}

void addToDB(int N, Database* db) {
	//updates a database with new states

	//set parameters
  double DT = 0.01;     //discretize samples into time chunks DT
	int pot = POTENTIAL; 
	int rho = RANGE;
	double beta = BETA;
	int method = 1;       //solve sde with em
	int steps = 200;      //number of time steps to take

	//initialize interaction matrices
	int* types = new int[N];
	int numTypes = readDesignFile(N, types);
	int numInteractions = numTypes*(numTypes+1)/2;
	double* kappa = new double[numInteractions];
	bd::readKappaFile(numInteractions, kappa);
	std::map<std::pair<int,int>, double> kmap; 
	bd::makeKappaMap(numTypes, kappa, kmap);
	int* P = new int[N*N];
	double* E = new double[N*N];
	bd::fillP(N, types, P, E, kmap);

	//set the initial and final state storage
	double* X = new double[DIMENSION*N];
	bd::setupChain(X, N);
	for (int i = 0; i < DIMENSION*N; i++) printf("%f\n", X[i]);

	//set up an adjacency matrix
	int* AM = new int[N*N];
	for (int i = 0; i < N*N; i++) {
		AM[i] = 0;
	}

	//create vector of new states
	std::vector<State> new_states;
	int count = 1;

	//do the time evolution
	for (int i = 0; i < steps; i++) {
		//solve sde
		solveSDE(X, N, DT, rho, beta, E, P, method, pot);

		//check if the state changed from previous step
		//get adjacency matrix for current state
		getAdj(X, N, AM);
		//printAM(N, AM);

		//check if this state has been seen before
		int state = searchDB(N, db, new_states, AM);
		std::cout << state << "\n";
		if (state == -1) {
			//add the new state to vector
			addState(N, X, AM, new_states);
			printf("Found new state. Total found this run: %d\n", count);
			count++;
		}

		if (i % (steps / 100) == 0) {
			printf("Finished step %d of %d\n", i, steps);
		}
	}

	//construct a new database
	int num_states_old = db->getNumStates();
	int num_states_new = new_states.size();
	int total = num_states_old + num_states_new;

	Database* newDB = new Database(N, total);

	//copy the first states from the old db
	for (int i = 0; i < num_states_old; i++) {
		(*newDB)[i] = (*db)[i];
	}

	//copy the rest from the vector
	for (int i = 0; i < num_states_new; i++) {
		(*newDB)[i+num_states_old] = new_states[i];
	}

	//print the new db
	std::string out = "N" + std::to_string(N) + "DBupdate.txt";
	std::ofstream out_str(out);
	out_str << *newDB;

	//delete memory
	delete []X; delete []AM;
	delete []types;
	delete []P; delete []E; delete []kappa;
	delete newDB;

}

















/******************************************************************/
/**************** Functions to sample MFPT    *********************/
/******************************************************************/


double sampleSTD(double* X, int n) {
	//evaluate the sample standard deviation of data in X
	double mean = 0; double std = 0;
	for (int i = 0; i < n; i++) mean += X[i];
	mean /= n;
	for (int i = 0; i < n; i++) std += (X[i]-mean)*(X[i]-mean);
	std /= (n-1); 
	return sqrt(std);
}


bool findMatrix(int* M, int* old, int old_bonds, int N, Database* db, int& timer, int& reset, int& reflect, 
																					int& new_state) {
	//find a state by transition matrix in the database

	for (int i = 0; i < (*db).getNumStates(); i++) {
		extractAM(N, i, old, db);
		bool S = checkSame(old, M, N); 
		if (S) {//found matrix
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
	bool C = checkConnected(M, N);
	if (!C) {//not connected
		reset = 1; timer -= 1;
		delete []M;
		return;
	}

	//compare old and new state, check if same
	int* old = new int[N*N]; for (int i = 0; i < N*N; i++) old[i]=0;
	int old_bonds = (*db)[state].getBonds();
	extractAM(N, state, old, db);
	bool S = checkSame(old, M, N);

	if (S) {//same matrix
		delete []M; delete []old;
		return;
	}
	else {//not the same, find matrix in database
		S = findMatrix(M, old, old_bonds, N, db, timer, reset, reflect, new_state);
		if (!S) {
			//may output an unphysical state. refine with newton, check again
			refine(N, X, M);
			S = findMatrix(M, old, old_bonds, N, db, timer, reset, reflect, new_state);
			if (!S) { //state still not found after refine, ignore this sample
				reset = 1; timer = 0;
				printf("State not found in database\n\n\n\n\n");
			}
		}
		delete []M; delete []old;
		return;
	}
}

void updatePM(int new_state, std::vector<Pair>& PM) {
	//find pair with index = state and increment value by 1

	int i;

	for (i = 0; i < PM.size(); i++) {
		if (PM[i].index == new_state) {
			PM[i].value += 1;
			break;
		}
	}

	if (i == PM.size()) {//this state is being hit for the first time
		PM.push_back(Pair(new_state,1));
	}

}

void runTrajectoryChain(double* X, int pot, Database* db, int state, int samples, int N, 
	double DT, int rho, double* E, double beta, int* P, int method, int& Num, 
																											int& Den, std::vector<Pair>& PM ) {
	//run the trajectory, update mfpt estimates

	//intiailize temp storage and set parameters
	double* temp = new double[2*N]; memcpy(temp, X, 2*N*sizeof(double));
	int reset; int reflect; int new_state = state; int hit = 0; int max_it = 100;
	int timer = 0; pot = 0;

	//solve sde and update
	for (int i = 0; i < max_it; i++) {
		reset = 0; reflect = 0;
		//solve SDE
		solveSDE(X, N, DT, rho, beta, E, P, method, pot);
		//check if state changed
		checkState(X, N, state, new_state, db, timer, reset, reflect);
		if (reflect == 0 && reset == 0) {//no hit, proceed
			memcpy(temp, X, 2*N*sizeof(double)); // copy x to temp
		}
		else if (reflect == 1) {//hit new state, update estimates
			Den += timer; Num += timer*(timer+1)/2.0; 
			//PM[new_state] += 1;
			updatePM(new_state, PM);
			memcpy(X, temp, 2*N*sizeof(double));//copy temp to x -> reset step
			timer = 0; hit +=1;
			//printf("New state %d", new_state);
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


void runTrajectoryMFPT(double* X, int pot, Database* db, int state, int samples, int N, 
	double DT, int rho, double* E, double beta, int* P, int method, int& Num, 
																											int& Den, std::vector<Pair>& PM ) {
	//run the trajectory, update mfpt estimates

	//intiailize temp storage and set parameters
	double* temp = new double[2*N]; memcpy(temp, X, 2*N*sizeof(double));
	int reset; int reflect; int new_state = state; int hit = 0; int max_it = 10*samples;
	int timer = 0; 

	//solve sde and update
	for (int i = 0; i < max_it; i++) {
		reset = 0; reflect = 0;
		//solve SDE
		solveSDE(X, N, DT, rho, beta, E, P, method, pot);
		//check if state changed
		checkState(X, N, state, new_state, db, timer, reset, reflect);
		if (reflect == 0 && reset == 0) {//no hit, proceed
			memcpy(temp, X, 2*N*sizeof(double)); // copy x to temp
		}
		else if (reflect == 1) {//hit new state, update estimates
			Den += timer; Num += timer*(timer+1)/2.0; 
			//PM[new_state] += 1;
			updatePM(new_state, PM);
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


void equilibrate(double* X, int pot, Database* db, int state, int eq, int N, double DT, 
													int rho, double* E, double beta, int* P, int method) {
	//perform eq steps to equilibrate the trajectory. do not record data

	//initialize temp storage and set parameters;
	double* temp = new double[2*N]; memcpy(temp, X, 2*N*sizeof(double));
	int reset; int reflect; int new_state = state; int timer = 0;

	//solve the SDE and update
	for (int i = 0; i < eq; i++) {
		reset = 0; reflect = 0;
		//solve the sde
		solveSDE(X, N, DT, rho, beta, E, P, method, pot);
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
	int pot = 1;  //set potential. 0 = morse, 1 = LJ
	int method = 1; //solve SDEs with EM
	int samples = 2000; //number of hits per walker for estimator
	int eq = 200; //number of steps to equilibrate for

	//quantities to update - old estimates
	int num_states = db->getNumStates();
	int num = (*db)[state].getNumerator(); 
	int den = (*db)[state].getDenominator(); 
	std::vector<Pair> pm = (*db)[state].getP();

	//quantities to update - new estimates
	int NUM = 0; int DEN = 0;
	std::vector<Pair> PM; std::vector<Pair> PMshare;
	double mfpt = 0;
	double sigma = 0;

	//output start message
	printf("Beginning MFPT Estimator for state %d out of %d.\n", state, num_states);

	//setup simulation
	double Eh = stickyNewton(8, rho, Kh, beta); //get energy corresponding to kappa
	//initialize interaction matrices
	int* P = new int[N*N]; double* E = new double[N*N];
	setupSimMFPT(N, Eh, P, E);

	//store mfpt estimates on each thread to get standard deviation
	double* mfptSamples; int num_threads;

	for (int i = 0; i < PM.size(); i++) printf("0 thread has:\n %d, %f\n", PM[i].index, PM[i].value);

	//open parallel region
	//#pragma omp parallel reduction(+:NUM, DEN, PM[:num_states])
	#pragma omp parallel private(PM) shared(PMshare) reduction(+:NUM, DEN)
	{

	//initialize samples
	num_threads = omp_get_num_threads();
	mfptSamples = new double[num_threads];

	//get starting structures
	const Cluster& c = (*db)[state].getRandomIC();

	//cluster structs to arrays
	double* X = new double[DIMENSION*N];
#if (DIMENSION == 2) 
	c.makeArray2d(X, N);
#endif
#if (DIMENSION == 3)
	c.makeArray3d(X, N);
#endif

	//equilibrate the trajectories
	equilibrate(X, pot, db, state, eq, N, DT, rho, E, beta, P, method);

	//run BD
	runTrajectoryMFPT(X, pot, db, state, samples, N, DT, rho, E, beta, P, method, NUM, DEN, PM );

	//store samples
	if (DEN != 0) {
		mfptSamples[omp_get_thread_num()] = (NUM * DT) / DEN;
	}
	else {
		mfptSamples[omp_get_thread_num()] = 0;
	}

	if (omp_get_thread_num() == 0) {
		PMshare = PM;
	}
	//do update on PM vectors - need barrier
	#pragma omp barrier
	#pragma omp critical
	{
		if (omp_get_thread_num() != 0) {
			combinePairs(PMshare, PM);
		}
	}

	//free cluster memory
	delete []X;

	//end parallel region
	}

	//combine estimates - if any samples were found (11->12 state)
	if (DEN != 0) {
		//combine estimates
		num += NUM; den += DEN;
		combinePairs(PMshare, pm); //PMshare has the updated info
		mfpt = (num * DT) / den;
		sigma = sampleSTD(mfptSamples, num_threads);

		//make a Z vector with same num of elements as P
		std::vector<Pair> Z; 
		for (int i = 0; i < PMshare.size(); i++) {
			Z.push_back(Pair(PMshare[i].index, 0));
		}

		//update database
		(*db)[state].mfpt = mfpt;
		(*db)[state].num = num;
		(*db)[state].denom = den;
		(*db)[state].num_neighbors = PMshare.size();
		(*db)[state].P = PMshare;
		(*db)[state].Z = Z;
		(*db)[state].Zerr = Z;
		(*db)[state].sigma = sigma;
	}

	//print out final estimates - debug
	/*
	printf("%Numerator = %d, Denominator = %d, \n", num, den);
	double sum = 0;
	for (int i = 0; i < PMshare.size(); i++) {
		printf("State = %d, visits = %f\n", PMshare[i].index, PMshare[i].value);
		sum +=PMshare[i].value;
	}
	printf("sum of hits = %f\n", sum);
	printf("Total Estimate = %f +- %f\n", mfpt, sigma);
	for (int i = 0; i < num_threads; i++) printf("MFPT estimate %d = %f\n", i, mfptSamples[i]);
	*/


	//free memory
	delete []E; delete []P; delete []mfptSamples;
	//delete []PM; delete []pm; 
}

void estimateChain(int N, int state, Database* db) {
	/*estimate mean first passage time starting in 
	chain state (must be given in state)*/


	//set parameters
	int rho = 40; double beta = 1; double DT = 0.01; int Kh = 1850;
	int pot = 1;  //set potential. 0 = morse, 1 = LJ
	int method = 1; //solve SDEs with EM
	int samples = 16667; //number of hits per walker for estimator

	//quantities to update - old estimates
	int num_states = db->getNumStates();
	int num = (*db)[state].getNumerator(); 
	int den = (*db)[state].getDenominator(); 
	std::vector<Pair> pm = (*db)[state].getP();

	//quantities to update - new estimates
	int NUM = 0; int DEN = 0;
	std::vector<Pair> PM; std::vector<Pair> PMshare;
	double mfpt = 0;
	double sigma = 0;

	//output start message
	printf("Beginning MFPT Estimator for linear chain.\n");

	//setup simulation
	double Eh = stickyNewton(8, rho, Kh, beta); //get energy corresponding to kappa
	//initialize interaction matrices
	int* P = new int[N*N]; double* E = new double[N*N];
	setupSimMFPT(N, Eh, P, E);
	printf("E = %f\n", Eh);
	

	//store mfpt estimates on each thread to get standard deviation
	double* mfptSamples; int num_threads;

	for (int i = 0; i < PM.size(); i++) printf("0 thread has:\n %d, %f\n", PM[i].index, PM[i].value);

	//open parallel region
	//#pragma omp parallel reduction(+:NUM, DEN, PM[:num_states])
	#pragma omp parallel private(PM) shared(PMshare) reduction(+:NUM, DEN)
	{

	//initialize samples
	num_threads = omp_get_num_threads();
	mfptSamples = new double[num_threads];

	double* X = new double[2*N];
	int mult = 1; int progress = 2000;



	//run BD
	for (int times = 0; times < samples; times++) {
		setupChain(X,N); 
		runTrajectoryChain(X, pot, db, state, 1, N, DT, rho, E, beta, P, method, NUM, DEN, PM );
		if (times % progress == 0) {
			printf("Thread %d generated sample %d.\n", omp_get_thread_num(), progress*mult);
			mult++;
		}
	}

	//store samples
	if (DEN != 0) {
		mfptSamples[omp_get_thread_num()] = DT * DEN / samples;
	}
	else {
		mfptSamples[omp_get_thread_num()] = 0;
	}

	if (omp_get_thread_num() == 0) {
		PMshare = PM;
	}
	//do update on PM vectors - need barrier
	#pragma omp barrier
	#pragma omp critical
	{
		if (omp_get_thread_num() != 0) {
			combinePairs(PMshare, PM);
		}
	}

	//free cluster memory
	delete []X;

	//end parallel region
	}

	//combine estimates - if any samples were found (11->12 state)
	if (DEN != 0) {
		//combine estimates
		num += NUM; den += DEN;
		combinePairs(PMshare, pm); //PMshare has the updated info
		mfpt = DT * den / (num_threads*samples);
		sigma = sampleSTD(mfptSamples, num_threads);

		//make a Z vector with same num of elements as P
		std::vector<Pair> Z; 
		for (int i = 0; i < PMshare.size(); i++) {
			Z.push_back(Pair(PMshare[i].index, 0));
		}

		//update database
		(*db)[state].mfpt = mfpt;
		(*db)[state].num = num;
		(*db)[state].denom = den;
		(*db)[state].num_neighbors = PMshare.size();
		(*db)[state].P = PMshare;
		(*db)[state].Z = Z;
		(*db)[state].Zerr = Z;
		(*db)[state].sigma = sigma;
	}

	//print out final estimates - debug
	/*
	printf("%Numerator = %d, Denominator = %d, \n", num, den);
	double sum = 0;
	for (int i = 0; i < PMshare.size(); i++) {
		printf("State = %d, visits = %f\n", PMshare[i].index, PMshare[i].value);
		sum +=PMshare[i].value;
	}
	printf("sum of hits = %f\n", sum);
	printf("Total Estimate = %f +- %f\n", mfpt, sigma);
	for (int i = 0; i < num_threads; i++) printf("MFPT estimate %d = %f\n", i, mfptSamples[i]);
	*/

	//free memory
	delete []E; delete []P; delete []mfptSamples;
	//delete []PM; delete []pm; 
}






/******************************************************************/
/**************** Functions to sample exit times ******************/
/******************************************************************/







double gyrationRadius(int N, double* X) {
	//compute the radius of gyration

	double* x = new double[N];
	double* y = new double[N];
	double* r2 = new double[N];

	double xCM = 0; double yCM = 0;

	for (int i = 0; i < N; i++) {
		x[i] = X[DIMENSION*i]; xCM += x[i];
		y[i] = X[DIMENSION*i+1]; yCM += y[i];
	}

	double rg = 0;
	for (int i = 0; i < N; i++) {
		r2[i] = (x[i]-xCM/N)*(x[i]-xCM/N) + (y[i]-yCM/N)*(y[i]-yCM/N);
		rg += r2[i];
	}

	rg /= N;

	//free memory
	delete []x; delete []y; delete []r2;

	return sqrt(rg);
}

double boop2d(int N, double* X) {
	//evaluate the 2-d bond orientational order parameter. (Nelson & Halperin)
	// phi_6 = 1/N sum_m 1/(NN_m) sum_NN_m e^(6i theta)
	//returns |phi_6|^2

	//construct the adjacency matrix
	int* M = new int[N*N]; for (int i = 0; i < N*N; i++) M[i]=0;
	getAdj(X, N, M);

	//init the boop - zero-d complex number
	std::complex<double> boop(0.0,0.0);

	//define i
	const std::complex<double> i(0.0,1.0);

	//loop over each particle
	for (int p1 = 0; p1 < N; p1++) {
		std::complex<double> summand = 0.0 + 0.0*i;
		int num_bonds = 0;
		//loop over each other particle, check for bonds, compute contribution
		for (int p2 = 0; p2 < N; p2++) {
			if (M[toIndex(p1,p2,N)] == 1) {
				//increment number of bonds
				num_bonds++;
				//compute vector lengths and dot product
				double x1 = X[DIMENSION*p1]; double y1 = X[DIMENSION*p1+1];
				double x2 = X[DIMENSION*p2]-x1; double y2 = X[DIMENSION*p2+1]-y1;;
				//double m1 = sqrt(x1*x1 + y1*y1); double m2 = sqrt(x2*x2 + y2*y2);
				//double dp = (x1*x2 + y1*y2) / (m1*m2);
				//double theta = acos(dp);
				double theta = atan(y2/x2);
				summand += std::exp(6.0 * i * theta);
			}
		}
		//divide by the number of bonds, add to boop
		summand /= num_bonds;
		boop += summand;
	}

	//divide boop by number of particles
	boop /= N;

	//compute the squared magnitude of the complex number boop
	double boop_norm2 = std::norm(boop);

	//free memory
	delete []M;

	//return the boop squared magnitude
	return boop_norm2;
}

double end2end(int N, double* X) {
	//end to end distance of a chain

	int p1 = 0; int p2 = N-1;
	double x1 = X[DIMENSION*p1]; double y1 = X[DIMENSION*p1+1];
	double x2 = X[DIMENSION*p2]; double y2 = X[DIMENSION*p2+1];
	double xD = x2-x1; double yD = y2-y1;

	return sqrt(xD*xD + yD*yD);

}

double rsa(int N, double* X) {
	//evaluate the relative shape anisotropy.
	//compute the gyration tensor, get eigenvalues, compute k^2

	double Sxx = 0; double Syy = 0; double Sxy = 0; 
	for (int p1 = 0; p1 < N; p1++) {
		for (int p2 = 0; p2 < N; p2++) {
			Sxx += (X[DIMENSION*p1] - X[DIMENSION*p2]) * (X[DIMENSION*p1] - X[DIMENSION*p2]);
			Syy += (X[DIMENSION*p1+1] - X[DIMENSION*p2+1]) * (X[DIMENSION*p1+1] - X[DIMENSION*p2+1]);
			Sxy += (X[DIMENSION*p1] - X[DIMENSION*p2]) * (X[DIMENSION*p1+1] - X[DIMENSION*p2+1]);
		}
	}

	Sxx /= (2.0*N); Sxy /= (2.0*N); Syy /= (2.0*N); 

	//compute the eigenvalues of the 2x2 symmetric matrix
	double D = (Sxx-Syy)*(Sxx-Syy) + 4*Sxy*Sxy;
	double l12 = 0.5 * (Sxx + Syy - sqrt(D));
	double l22 = 0.5 * (Sxx + Syy + sqrt(D));

	//compute the shape parameters
	double b = -0.5 * (l12 + l22);
	double c = l22 - l12;
	double Rg2 = l12 + l22;

	//compute the ksa
	double k2 = (b*b + 0.75*c*c) / (Rg2*Rg2);
	return k2;
}

void setupChain(double* X, int N, int which) {
	//initialize coordinates of a chain via some choice

	if (which == 0) {
		for (int i = 0; i < N; i++) {
			X[2*i] = i; X[2*i+1] = 0;
		}
	}
	else if (which == 1) { //1 bond
		for (int i = 0; i < N; i++) {
			X[2*i] = i; X[2*i+1] = 0;
		}
		X[10] = 3.5; X[11] = sqrt(3.0) / 2.0;
	}
}

void sampleFirstExit(int N, int initial, Database* db) {
	/*get samples of some quantity at the firste exit time, starting from a linear chain */

	//set parameters
	int rho = 40; double beta = 1; double DT = 0.01; int Kh = 1850;
	int pot = 0;  //set potential. 0 = morse, 1 = LJ
	int method = 1; //solve SDEs with EM
	int samples = 500; //number of samples to get

	//cutoff for qsd
	int t_cut = 30;

	//setup simulation
	double Eh = stickyNewton(8, rho, Kh, beta); //get energy corresponding to kappa
	//initialize interaction matrices
	int* P = new int[N*N]; double* E = new double[N*N];
	setupSimMFPT(N, Eh, P, E);
	printf("E = %f\n", Eh);

	//make a vector to store samples
	std::vector<double> q_samples;

	//run BD
	#pragma omp parallel 
	{
	//setup position storage
	double* X = new double[DIMENSION*N];
	double* temp = new double[DIMENSION*N];
	#pragma omp parallel for
	for (int times = 0; times < samples; times++) {

		printf("Running estimate %d\n", times+1);
		setupChain(X,N,0); 
		int state = initial;

		//intiailize temp storage and set parameters
	  memcpy(temp, X, N*DIMENSION*sizeof(double));
		int reset; int reflect; int new_state = state; int max_it = 10*samples;
		int timer = 0; 

		//these lines are debug to find the index of a state defined by coordinates
		//checkState(X, N, state, new_state, db, timer, reset, reflect);
		//printf("First state is %d\n", new_state);
		//abort();

		//solve sde and update
		for (int i = 0; i < max_it; i++) {
			reset = 0; reflect = 0;
			//solve SDE
			solveSDE(X, N, DT, rho, beta, E, P, method, pot);

			//check if state changed
			checkState(X, N, state, new_state, db, timer, reset, reflect);

			if (reflect == 0 && reset == 0) { //no hit, proceed
				memcpy(temp, X, DIMENSION*N*sizeof(double)); // copy x to temp

			}
			else if (reflect == 1) {//hit new state, get sample of quantity
				if (i > t_cut) {
					//double q = gyrationRadius(N, X);
					//double q = boop2d(N, X);
					double q = end2end(N, X);
					q_samples.push_back(q);
					std::cout << i << "\n";
					break;
				}
				else {
					memcpy(X, temp, DIMENSION*N*sizeof(double));//copy temp to x -> reset step
				}

			}
			else {//chain broke, reset previous config
				memcpy(X, temp, DIMENSION*N*sizeof(double));//copy temp to x -> reset step

			}
		}
	}
	delete []X; delete []temp;
	}

	//output the samples to a file
	std::ofstream ofile;
	ofile.open("fhtBD.txt");
	for (int i = 0; i < q_samples.size(); i++) {
		ofile << q_samples[i] << "\n";
	}
	ofile.close();
	
	//free memory
	delete []E; delete []P; 
}

void sampleSecondExit(int N, int initial, Database* db) {
	/*get samples of some quantity at the firste exit time, starting from a linear chain */

	//set parameters
	int rho = 40; double beta = 1; double DT = 0.01; int Kh = 1850;
	int pot = 0;  //set potential. 0 = morse, 1 = LJ
	int method = 1; //solve SDEs with EM
	int samples = 300; //number of samples to get

	//cutoff for qsd
	int t_cut = 500;

	//setup simulation
	double Eh = stickyNewton(8, rho, Kh, beta); //get energy corresponding to kappa
	//initialize interaction matrices
	int* P = new int[N*N]; double* E = new double[N*N];
	setupSimMFPT(N, Eh, P, E);
	printf("E = %f\n", Eh);

	//make a vector to store samples
	std::vector<double> q_samples;

	//run BD
	#pragma omp parallel
	{
	//setup position storage
	double* X = new double[DIMENSION*N];
	double* temp = new double[DIMENSION*N];
	for (int times = 0; times < samples; times++) {

		printf("Running estimate %d\n", times+1);
		setupChain(X,N); 
		int state = initial;

		//intiailize temp storage and set parameters
	  memcpy(temp, X, N*DIMENSION*sizeof(double));
		int reset; int reflect; int new_state = state; int max_it = 2000;
		int timer = 0; 

		//solve sde and update
		for (int i = 0; i < max_it; i++) {
			reset = 0; reflect = 0;
			//solve SDE
			solveSDE(X, N, DT, rho, beta, E, P, method, pot);

			//check if state changed
			checkState(X, N, state, new_state, db, timer, reset, reflect);
			int new_bonds = (*db)[new_state].getBonds();

			if (new_bonds < 7 && reset == 0) { //no hit, proceed
				memcpy(temp, X, DIMENSION*N*sizeof(double)); // copy x to temp
				state = new_state;
			}
			else if (new_bonds == 7 && reflect == 1) {//hit new state, get sample of quantity
				if (i > t_cut) {
					//double q = gyrationRadius(N, X);
					//double q = boop2d(N, X);
					double q = end2end(N, X);
					q_samples.push_back(q);
					std::cout << i << "\n";
					break;
				}
				else {
					memcpy(X, temp, DIMENSION*N*sizeof(double));//copy temp to x -> reset step
				}

			}
			else {//chain broke, reset previous config
				memcpy(X, temp, DIMENSION*N*sizeof(double));//copy temp to x -> reset step

			}
		}
	}
	delete []X; delete []temp;
	}

	//output the samples to a file
	std::ofstream ofile;
	ofile.open("fhtBD.txt");
	for (int i = 0; i < q_samples.size(); i++) {
		ofile << q_samples[i] << "\n";
	}
	ofile.close();
	
	//free memory
	delete []E; delete []P; 
}

void fillHydroData(int N, std::string filename, std::vector<std::vector<double>>& ics) {
	//take hydro data from file and put in vector

	std::ifstream in_str(filename);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		return;
	}

	double x0; //check if there is another line
	double x;
	std::vector<double> cluster;

	//fill the vector
	while (in_str >> x0) {
		cluster.clear();
		cluster.push_back(x0);
		for (int i = 0; i < N*DIMENSION-1; i++) {
			in_str >> x; cluster.push_back(x);
		}
		ics.push_back(cluster);
	}

}

void sampleSecondExit(int N, Database* db) {
	/*get samples of some quantity at the second exit time. Uses the hydrodynamics 
	  data at the first exit time as the initial condition.  */

	//set parameters
	int rho = 40; double beta = 1; double DT = 0.01; int Kh = 1850;
	int pot = 0;  //set potential. 0 = morse, 1 = LJ
	int method = 1; //solve SDEs with EM

	//get HD data
	int hydro = 0;   //set to 0 or 1 for ex/including short range hydrodynamics
	std::string base = "input/hydro/"; 
	//check if hydrodynamics were on or off, set parameters accordingly
	if (hydro == 0) { //hydro is off
	  base += "noHD/noHD1bond.txt";
	}
	else { //hydro is on
		base += "HD/HD1bond.txt";
	}
	std::vector<std::vector<double>> ics;
	fillHydroData(N, base, ics);
	int num_samples = ics.size();


	//cutoff for qsd
	int t_cut = 50;

	//setup simulation
	double Eh = stickyNewton(8, rho, Kh, beta); //get energy corresponding to kappa
	//initialize interaction matrices
	int* P = new int[N*N]; double* E = new double[N*N];
	setupSimMFPT(N, Eh, P, E);
	printf("E = %f\n", Eh);

	//make a vector to store samples
	//std::vector<double> q_samples;
	double* q_samples = new double[num_samples];
	for (int i = 0; i < num_samples; i++) {
		q_samples[i] = 0.0;
	}

	//run BD
	#pragma omp parallel
	{
	//setup position storage
	double* X = new double[DIMENSION*N];
	double* temp = new double[DIMENSION*N];
	#pragma omp for schedule(auto)
	for (int sample = 0; sample < num_samples; sample++) {

		std::vector<double> coordinates = ics[sample];
		for (int c = 0; c < N*DIMENSION; c++) {
			X[c] = coordinates[c];
		}
		int dummy = 0; int state;
		checkState(X, N, 1, state, db, dummy, dummy, dummy);
		printf("Sample %d on thread %d is starting in state %d\n", sample, omp_get_thread_num(), state);
		

		//intiailize temp storage and set parameters
	  memcpy(temp, X, N*DIMENSION*sizeof(double));
		int reset; int reflect; int new_state = state; int max_it = 2000;
		int timer = 0; 

		//solve sde and update
		for (int i = 0; i < max_it; i++) {
			reset = 0; reflect = 0;
			//solve SDE
			solveSDE(X, N, DT, rho, beta, E, P, method, pot);

			//check if state changed
			checkState(X, N, state, new_state, db, timer, reset, reflect);
			int new_bonds = (*db)[new_state].getBonds();

			if (new_bonds < 7 && reset == 0) { //no hit, proceed
				memcpy(temp, X, DIMENSION*N*sizeof(double)); // copy x to temp
				state = new_state;
			}
			else if (new_bonds == 7 && reflect == 1) {//hit new state, get sample of quantity
				if (i > t_cut) {
					//double q = gyrationRadius(N, X);
					//double q = boop2d(N, X);
					double q = end2end(N, X);
					//q_samples.push_back(q);
					q_samples[sample] = q;
					std::cout << i << "\n";
					break;
				}
				else {
					memcpy(X, temp, DIMENSION*N*sizeof(double));//copy temp to x -> reset step
				}

			}
			else {//chain broke, reset previous config
				memcpy(X, temp, DIMENSION*N*sizeof(double));//copy temp to x -> reset step

			}
		}
	}
	delete []X; delete []temp;
	}

	//output the samples to a file
	std::ofstream ofile;
	ofile.open("fhtBD.txt");
	for (int i = 0; i < num_samples; i++) {
		if (q_samples[i] != 0) {
			ofile << q_samples[i] << "\n";
		}
	}
	ofile.close();
	
	//free memory
	delete []E; delete []P; 
	delete []q_samples;
}







}