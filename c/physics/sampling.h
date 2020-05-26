#pragma once
#include "nauty.h"
#include "adjacency.h"
#include "pair.h"
#include <vector>
#include <random>
#include <chrono>
#include <eigen3/Eigen/Dense>

class RandomNo{
	std::mt19937 generator; 
	std::uniform_real_distribution<double> uDist;
	std::normal_distribution<double> gDist;

	public:
		RandomNo(unsigned int seed = std::random_device{}())
        : generator{seed},uDist{0,1},gDist{0,1} {}

    double getU() {
    	return uDist(generator);
    }
    double getG() {
    	return gDist(generator);
    }


};

namespace bd { 
class Database;

//brownian dynamics sampling section

//building a database of all states by sampling


//doing mfpt estimation with completed database
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
double sampleSTD(double* X, int n);
bool findMatrix(int* M, int* old, int old_bonds, int N, Database* db, int& timer, 
	int& reset, int& reflect, int& new_state);


//sampling quantities at exit times
void sampleFirstExit(int N, int state, Database* db);
void sampleSecondExit(int N, int state, Database* db);
void sampleSecondExit(int N, Database* db);

//order parameters
double gyrationRadius(int N, double* X);
double boop2d(int N, double* X);
double end2end(int N, double* X);
double rsa(int N, double* X);



}

namespace mcm {
class Database;

//monte carlo on manifolds sampling section

//trimer test functions
void runTestTrimer(int N, double* X, int* M, int num_samples);
double getBondAngle(Eigen::VectorXd x);

//mcmc algorithm functions
bool getSample(int N, int* M, int df, int b, int d, Eigen::VectorXd& x, RandomNo* rng);
bool project(int N, int* M, int b, Eigen::VectorXd z, Eigen::MatrixXd Qz, Eigen::VectorXd& a);
bool checkInequality(int N, int* M, Eigen::VectorXd y);
double getResidual(int N, int* M, int b, Eigen::VectorXd p);
void QRortho(Eigen::MatrixXd& Qx, Eigen::MatrixXd& Q, int d);
int getNumBonds(int N, int* M);
double evalDensity(Eigen::VectorXd v, int d);
void proposeTan(Eigen::MatrixXd Q, Eigen::VectorXd& v, int d, RandomNo* rngee);
void evalConstraint(int N, int* M, Eigen::VectorXd& q, Eigen::VectorXd x);
void makeJ(int N, int* M, Eigen::MatrixXd& J, Eigen::VectorXd x);
void makeQx(int N, int* M, Eigen::MatrixXd& Qx, Eigen::VectorXd x);
double getMH(Eigen::MatrixXd Qx, Eigen::MatrixXd Qy,Eigen::VectorXd x, 
						 Eigen::VectorXd y, Eigen::VectorXd v, Eigen::VectorXd vr, int d, int b);
double fEval(int b, Eigen::MatrixXd Q);
double getJacobian(Eigen::MatrixXd Q);

//chain mfpt estimator functions
void sampleStats(std::vector<double> X, double& M, double& V);
void sampleStats(double* X, int N, double& M, double& V);
void minVarEstimate(int sampleSize, double* means, double* variances, double& M, double& V);
double getTime(int N, Eigen::VectorXd x, Eigen::VectorXd x0);
bool findMatrix(int* M, int* old, int old_bonds, int N, bd::Database* db, 
								bool& reset, int& new_state);
void checkState(int N, Eigen::VectorXd x, int state, bd::Database* db,
							  bool& reset, int& new_state);
void equilibrate(double* X, bd::Database* db, int state, int N, int* M, RandomNo* rngee);
double getSampleMFPT(double* X, bd::Database* db, int state, int N, int* M,
	std::vector<bd::Pair>& PM, RandomNo* rngee);
void getSamplesMFPT(double* X, bd::Database* db, int state, int N, int* M,
	std::vector<bd::Pair>& PM, std::vector<double>& mfptVec, RandomNo* rngee);
void estimateMFPTreset(int N, int state, bd::Database* db, RandomNo* rngee);
void estimateMFPTreflect(int N, int state, bd::Database* db, RandomNo* rngee);




}


