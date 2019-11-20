#pragma once
#include <vector>
#include <map>
#include <deque>
namespace bd { 
class Database; 

//general functions
void getBondTypes(int N, int* particleTypes, Database* db, std::vector<int> targets);
int setTypes(int N, int* particleTypes, int IC);
double getStickyProduct(int N, int state,  Database* db, int* particleTypes, 
												std::map<std::pair<int,int>,double> kappa);
void reweight(int N, int num_states, Database* db, int* particleTypes, double* eq,
						  std::map<std::pair<int,int>,double> kappa);
void initKappaVals(int numInteractions, double* kappaVals);
void constructSurfaceTOY(int N, Database* db, int initial, int target, bool useFile);
void allPerms(int N, std::deque<std::string>&);
void distinctPerms(int N, std::deque<std::string>&);
void checkPerm(std::string , std::deque<std::string>&);
void checkPositive(int numInteractions, double* kappaVals); 
void applyMax(int numInteractions, double* kappaVals, double M);

//hitting probability optimization
double getHitProbability(int num_states, int initial, std::vector<int> targets, double* U);
double solveAbsorbProblem(int initial, double* kappaVals, Database* db, int* particleTypes, 
								 					double* Tconst, std::vector<int> ground, std::vector<int> targets);

double computeGradHP(int initial, int numInteractions, double* kappaVals, Database* db, 
								 int* particleTypes, double* Tconst, std::vector<int> ground, 
								 std::vector<int> targets, double* g);

void hittingProbMaxTOY(int N, Database* db, int initial, int target, bool useFile);
void hittingProbMaxTOYperms(int N, Database* db, int initial, int target);

//equilibrium probability optimization
void eqProbMaxTOY(int N, Database* db, int initial, int target, bool useFile);
void eqProbMaxTOYperms(int N, Database* db, int initial, int target);
double computeGradEQ(int initial, int numInteractions, double* kappaVals, Database* db, 
								 int* particleTypes, std::vector<int> targets, double* g);
double getEqProb(int initial, double* kappaVals, Database* db, 
								 int* particleTypes, std::vector<int> targets);

//average transition rate optimization
void rateMaxTOY(int N, Database* db, int initial, int target, bool useFile);
void rateMaxTOYperms(int N, Database* db, int initial, int target);
double computeGradRate(int initial, int numInteractions, double* kappaVals, Database* db, 
								 int* particleTypes, double* Tconst, std::vector<int> targets, double* g);
double getRate(int initial, double* kappaVals, Database* db, int* particleTypes, 
							 double* Tconst, std::vector<int> targets);

//hitting probability optimization with rate constraints
void hittingProbMaxTOYc(int N, Database* db, int initial, int target, bool useFile);



//sampling functions
int findMatrix(int N, double* X, Database* db);
void estimateHittingProbability(int N, Database* db, int target);









}