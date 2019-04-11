#pragma once
#include "point.h"
#include <../Eigen/Dense>
#include "nauty.h"
#include "database.h"
#include <string>
#include <vector>

namespace bd {

class Database;

//adjacency matrix stuff
int toIndex(int r, int c, long m);
void index2ij(int index, int N, int& i, int& j);
int checkConnected(int* M, int N);
int checkSame(int* M1, int* M2, int N);
void getAdj(double* X, int N, int* M);
void findIsomorphic(int N, int num_states, int state, Database* db, std::vector<int>&);
void buildNautyGraph(int N, int M, int state, Database* db, graph* g);
bool checkIsomorphic(int N, int M, graph* g1, graph* g2);

//evaluate euclidean distance 
double euDist(double* particles, int i, int j, int N, double* Z);
//go from cluster array to particle array
void c2p(double* cluster, double* particles, int N);

}