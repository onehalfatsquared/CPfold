#pragma once
#include <vector>
#include <map>
namespace bd { 
class Database; 

//general functions
int setTypes(int N, int* particleTypes, int IC);
double getStickyProduct(int N, int state,  Database* db, int* particleTypes, 
												std::map<std::pair<int,int>,double> kappa);
void reweight(int N, int num_states, Database* db, int* particleTypes, double* eq,
						  std::map<std::pair<int,int>,double> kappa);
void initKappaVals(int numInteractions, double* kappaVals);
void constructSurfaceTOY(int N, Database* db, int initial, int target, bool useFile);
double getHitProbability(int num_states, int initial, std::vector<int> targets, double* U);

//optimization functions
double solveAbsorbProblem(int initial, double* kappaVals, Database* db, int* particleTypes, 
								 					double* Tconst, std::vector<int> ground, std::vector<int> targets);

double computeGradHP(int initial, int numInteractions, double* kappaVals, Database* db, 
								 int* particleTypes, double* Tconst, std::vector<int> ground, 
								 std::vector<int> targets, double* g);

void hittingProbMaxTOY(int N, Database* db, int initial, int target, bool useFile);
void hittingProbMaxTOYperms(int N, Database* db, int initial, int target);

//sampling functions
int findMatrix(int N, double* X, Database* db);
void estimateHittingProbability(int N, Database* db, int target);









}