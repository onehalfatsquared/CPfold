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
void equilibrate();
void getSampleMFPT();
void minVarEstimate(int sampleSize, double* means, double* variances);
void estimateMFPT(int N, int state, bd::Database* db);


}


