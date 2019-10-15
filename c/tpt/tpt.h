#pragma once
#include <vector>
namespace bd { 
class Database; 

//general functions
//perform tpt on initial and target states
void performTPT(int N, int initial, int target, Database* db, bool getIso);

//transition matrix functions
//fill in rate matrix with data from mfpt estimates
void createTransitionMatrix(double* T, int num_states, Database* db, 
														std::vector<int>& endStates);
//create the invariant measure for the current problem
void createMeasure(int num_states, Database* db, double* eq, double kappa);
//create configurational partition fn for each cluster
void computePartitionFn(int num_states, Database* db, double* Z); 
//edit rate matrix to satisfy detailed balance
void satisfyDB(double* T, int num_states, Database* db, double* eq);
//fill in diagonal of rate matrix
void fillDiag(double* T, int num_states);
//create probability transition matrix
void createProbabilityMatrix(double* T, int num_states, double* P); 

//tpt functions
//set up and solve dirchlet problem for committor function
void computeCommittor(double* q, double* T, int num_states, int initial, std::vector<int>);
//compute the probability flux from committor
void computeFlux(int num_states, double* q, double* T, double* eq, double* flux);
//compute free energy of each cluster - from configurational partition fn
void computeFreeEnergy(int num_states, double* Z, double* F); 
//compute hitting probabilities for end states in reversible case
void computeHittingProbability(double* P, int num_states, std::vector<int> endStates, 
															 double* U);

























}