#pragma once
namespace bd { 
class Database;

//particle and potential stuff
//evaluate euclidean distance 
double euDist(double* particles, int i, int j, int N, double* Z);
//go from cluster array to particle array
void c2p(double* cluster, double* particles, int N);
//evaluate the morse potential between 2 particles
double morseP(double r, double rho, double E);
//evaluate total potential energy of system of particles
double morseEval(double* particles, int rho, double* E, int N, int* P);
//evaluate gradient of morse potential
void morseGrad(double* particles, int rho, double* E, int N, int* P, double* g);
//evaluate the sticky parameter for morse potential
void stickyF(double E, double rho, double k0, double& f, double& fprime);
//use newtons method to find E from kappa
double stickyNewton(double E, double rho, double k0);

//adjacency matrix stuff
int toIndex(int r, int c, long m);
void index2ij(int index, int N, int& i, int& j);
int checkConnected(int* M, int N);
int checkSame(int* M1, int* M2, int N);
void getAdj(double* X, int N, int* M);

//time integration stuff
//use em scheme to solve sde
void EM(double* X0, int N, int Nt, double k, 
					int rho, double* E, int* P, double beta);
//solve sde system
void solveSDE(double* X0, int N, double T, int rho, double beta,
							double* E, int* P, int method);

//sampling stuff
void estimateMFPT(int N, int state, Database* db);
void setupSim(int N, double Eh, int*& P, double*& E);
void equilibrate(double* X, Database* DB, int state, int eq, int N, double DT,
													int rho, double* E, double beta, int* P, int method);
void runTrajectory(double* X, Database* DB, int state, int samples, int N, 
	double DT, int rho, double* E, double beta, int* P, int method, int& Num, 
																									int& Den, int* PM );
void checkState(double* X, int N, int state, int new_state, Database* db, int& timer,
							 int& reset, int& reflect);
void extractAM(int N, int state, int* AM, Database* db);



}