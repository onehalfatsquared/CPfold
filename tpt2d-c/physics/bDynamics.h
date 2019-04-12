#pragma once
#include "nauty.h"
#include "adjacency.h"
#include "pair.h"
#include <vector>
#include <../Eigen/Dense>
namespace bd { 
class Database;

//particle and potential stuff
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


//time integration stuff
//use em scheme to solve sde
void EM(double* X0, int N, int Nt, double k, 
					int rho, double* E, int* P, double beta);
//solve sde system
void solveSDE(double* X0, int N, double T, int rho, double beta,
							double* E, int* P, int method);

//sampling stuff
void estimateMFPT(int N, int state, Database* db);
void estimatePartitionFn(int N, int state, Database* db);
void setupSimMFPT(int N, double Eh, int*& P, double*& E);
void equilibrate(double* X, Database* DB, int state, int eq, int N, double DT,
													int rho, double* E, double beta, int* P, int method);
void runTrajectoryMFPT(double* X, Database* DB, int state, int samples, int N, 
	double DT, int rho, double* E, double beta, int* P, int method, int& Num, 
																									int& Den, std::vector<Pair>& PM );
void updatePM(int new_state, std::vector<Pair>& PM); 
void checkState(double* X, int N, int state, int& new_state, Database* db, int& timer,
							 int& reset, int& reflect);
void refine(int N, double* X, int* M);
void makeNM(int N, int* M, Eigen::VectorXd x, Eigen::MatrixXd& J, Eigen::VectorXd& F);
void extractAM(int N, int state, int* AM, Database* db);
double sampleSTD(double* X, int n);
bool findMatrix(int* M, int* old, int old_bonds, int N, Database* db, int& timer, 
	int& reset, int& reflect, int& new_state);



}