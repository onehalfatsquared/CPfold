#pragma once
#include "nauty.h"
#include "adjacency.h"
#include "pair.h"
#include <fstream>
#include <map>
#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
namespace bd { 
class Database;

//general stuff
void makeKappaMap(int numTypes, double* kappaVals, 
									std::map<std::pair<int,int>,double>& kappa);
void fillP(int N, int* particleTypes, int* P, double* E,
					 std::map<std::pair<int,int>,double>& kappa);
void readKappaFile(int numInteractions, double* kappa);
int readDesignFile(int N, int* particleTypes);

//particle and potential stuff
//evaluate the morse potential between 2 particles
double morseP(double r, double rho, double E);
//evaluate total potential energy of system of particles
double morseEval(double* particles, int rho, double* E, int N, int* P);
//evaluate gradient of morse potential
void morseGrad(double* particles, int rho, double* E, int N, int* P, double* g);
void morseGradR(double* particles, int rho, double* E, int N, int* P, double* g);
//evlauate lj potential
double ljP(double r, double rho, double E);
//evaluate total potential energy of system of particles
double ljEval(double* particles, int rho, double* E, int N, int* P);
//evaluate gradient of morse potential
void ljGrad(double* particles, int rho, double* E, int N, int* P, double* g);
//evaluate the sticky parameter for morse potential
void stickyF(double E, double rho, double beta, double k0, double& f, double& fprime);
//use newtons method to find E from kappa
double stickyNewton(double E, double rho, double k0, double beta);


//time integration stuff
//init the chain
void setupChain(double* X, int N);
//print the cluster
void printCluster(double* X, int N);
//use em scheme to solve sde
void EM(double* X0, int N, int Nt, double k, 
					int rho, double* E, int* P, double beta, int pot);
//solve sde system
void solveSDE(double* X0, int N, double T, int rho, double beta,
							double* E, int* P, int method, int pot);





}