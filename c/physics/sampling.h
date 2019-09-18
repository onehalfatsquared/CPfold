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

void runSampler(int N, double* X, int* M);

bool project(int N, int* M, int b, Eigen::VectorXd z, Eigen::MatrixXd Qz, Eigen::VectorXd& a);
bool checkInequality(int N, int* M, Eigen::VectorXd y);
double getResidual(int N, int* M, int b, Eigen::VectorXd p);
void QRortho(Eigen::MatrixXd& Qx, Eigen::MatrixXd& Q, int d);
int getNumBonds(int N, int* M);
double evalDensity(Eigen::VectorXd v, int d);
void proposeTan(Eigen::MatrixXd Q, Eigen::VectorXd& v);
void evalConstraint(int N, int* M, Eigen::VectorXd& q, Eigen::VectorXd x);
void makeJ(int N, int* M, Eigen::MatrixXd& J, Eigen::VectorXd x);
void makeQx(int N, int* M, Eigen::MatrixXd& Qx, Eigen::VectorXd x);
double getMH(Eigen::VectorXd x, Eigen::VectorXd y, Eigen::VectorXd v, Eigen::VectorXd vr, int d, int b);
double fEval(int b);


}

