#pragma once
#include "point.h"
#include <eigen3/Eigen/Dense>
#include "nauty.h"
#include "database.h"
#include <string>
#include <vector>

namespace bd {

class Database;

//adjacency matrix stuff
int toIndex(int r, int c, long m);
void index2ij(int index, int N, int& i, int& j);
bool checkConnected(int* M, int N);
bool checkSame(int* M1, int* M2, int N);
void getAdj(double* X, int N, int* M);
void getAdjCut(double* X, int N, int* M, double cut);
void findIsomorphic(int N, int num_states, int state, Database* db, std::vector<int>&);
void buildNautyGraph(int N, int M, int state, Database* db, graph* g);
bool checkIsomorphic(int N, int M, graph* g1, graph* g2);
void extractAM(int N, int state, int* AM, Database* db);
void printAM(int N, int* AM);
void refine(int N, double* X, int* M);
void makeNM(int N, int* M, Eigen::VectorXd x, Eigen::MatrixXd& J, Eigen::VectorXd& F);

//evaluate euclidean distance 
double euDist(double* particles, int i, int j, int N, double* Z);
//go from cluster array to particle array
void c2p(double* cluster, double* particles, int N);

}