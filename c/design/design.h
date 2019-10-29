#pragma once
#include <vector>
#include <map>
namespace bd { 
class Database; 


int setTypes(int N, int* particleTypes, int IC);
int readDesignFile(int N, int* particleTypes);
double getStickyProduct(int N, int state,  Database* db, int* particleTypes, 
												std::map<std::pair<int,int>,double> kappa);
void reweight(int N, int num_states, Database* db, int* particleTypes, double* eq,
						  std::map<std::pair<int,int>,double> kappa);
void initKappaVals(int numInteractions, double* kappaVals);
void makeKappaMap(int numTypes, double* kappaVals, 
									std::map<std::pair<int,int>,double>& kappa);
void constructSurfaceTOY(int N, Database* db, int initial, int target, bool useFile);










}