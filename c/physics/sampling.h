#pragma once
#include "nauty.h"
#include "adjacency.h"
#include "pair.h"
#include <vector>
#include <eigen3/Eigen/Dense>
namespace bd { 
class Database;

//brownian dynamics sampling stuff
void estimateMFPT(int N, int state, Database* db);
void estimateChain(int N, int state, Database* db);
void setupSimMFPT(int N, double Eh, int*& P, double*& E);
void equilibrate(double* X, int pot, Database* DB, int state, int eq, int N, double DT,
													int rho, double* E, double beta, int* P, int method);
void runTrajectoryMFPT(double* X, int pot, Database* DB, int state, int samples, int N, 
	double DT, int rho, double* E, double beta, int* P, int method, int& Num, 
																									int& Den, std::vector<Pair>& PM );
void runTrajectoryChain(double* X, int pot, Database* db, int state, int samples, int N, 
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

namespace mcm {

//monte carlo on manifolds sampling stuff





}

