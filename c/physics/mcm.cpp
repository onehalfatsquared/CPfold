#include <math.h>
#include <random>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <chrono>
#include "bDynamics.h"
#include "database.h"
#include "adjacency.h"
#include "../defines.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/QR>
#include <omp.h>
namespace mcm {

void project(Eigen::VectorXd z, Eigen::MatrixXd Qz) {
	//project the point onto the manifold using newtons method

}

double evalDensity(Eigen::VectorXd v, int d) {
	//evaluate the gaussian d dimensional probability desnity at v
	double V = v.squaredNorm();  //get the squared 2 norm of v
	double Z = 1;               //compute normalizing const with for loop
	for (int i = 0; i < d; i++) {
		Z /= (sqrt(2*M_PI)*SIG);
	}
	return Z*exp(-V/(2*SIG*SIG));
}

void proposeTan(Eigen::MatrixXd Q, Eigen::VectorXd& v, int d) {
	//generate a proposal move tangent to M - isotropic gaussian

	//initialize the random number generator - Normal(0,1) //fix seed
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<double> distribution(0.0,1.0);

	for (int i = 0; i < d; i++) {
		v += SIG*distribution(generator)*Q.col(i);
	}

}

void QRortho(Eigen::MatrixXd& Qx, Eigen::MatrixXd& Q, int d) {
	//get the last d cols of Q matrix from a QR decomp of Qx

	Eigen::MatrixXd Qfull = Qx.colPivHouseholderQr().matrixQ();
	Q = Qfull.rightCols(d);
}

void makeQx(int N, int* M, Eigen::MatrixXd& Qx, Eigen::VectorXd x) {
	//construct the (transpose of the) jacobian of the constraint matrix

	double XD, YD, ZD;
	int count = 0;

	//loop over and construct system
	for (int i = 0; i <N; i ++) {
		for (int j = i+1; j < N; j++) {
			if (M[bd::toIndex(i,j,N)] == 1) {
				//compute distances
				XD = x(DIMENSION*i) - x(DIMENSION*j);
				YD = x(DIMENSION*i+1) - x(DIMENSION*j+1);
#if (DIMENSION == 3) 
				ZD = x(DIMENSION*i+2) - x(DIMENSION*j+2);
#endif 

				Qx(count, DIMENSION*i) = 2*XD; Qx(count, DIMENSION*j) = -2*XD;
				Qx(count, DIMENSION*i+1) = 2*YD; Qx(count, DIMENSION*j+1) = -2*YD;
#if (DIMENSION == 3) 
				Qx(count, DIMENSION*i+2) = 2*ZD; Qx(count, DIMENSION*j+2) = -2*ZD;
#endif

				count +=1;
			}
		}
	}

	//get the transpose of Qx
	Qx.transposeInPlace();
}

void evalConstraint(int N, int* M, Eigen::VectorXd& q, Eigen::VectorXd x) {
	//evaluate the lhs constraint function

	double XD, YD;
	double ZD = 0;
	int count = 0;

	//loop over and construct system
	for (int i = 0; i <N; i ++) {
		for (int j = i+1; j < N; j++) {
			if (M[bd::toIndex(i,j,N)] == 1) {
				//compute distances
				XD = x(DIMENSION*i) - x(DIMENSION*j);
				YD = x(DIMENSION*i+1) - x(DIMENSION*j+1);
#if (DIMENSION == 3) 
				ZD = x(DIMENSION*i+2) - x(DIMENSION*j+2);
#endif 

				//evaluate the distance constraint
				q(count) = XD*XD + YD*YD + ZD*ZD -1;

				count +=1;
			}
		}
	}
}

int getNumBonds(int N, int* M) {
	//return the number of bonds from an adjacency matrix

	int b = 0;
	for (int i = 0; i < N*N; i++) {
		b += M[i];
	}
	b /= 2;
	return b;
}

void setup(int N, double* X, int* M) {
	//do stuff //todo

	//get the number of bonds in the configuration, get DoF
	int df = DIMENSION*N;           //num degrees of freedom in ambient system
	int b = getNumBonds(N,M);      //num bonds
	int d = df-b;                  //effective dimension of system

	//setup the matrices/vectors with eigen
	Eigen::MatrixXd Qx(b,df); 
	Eigen::VectorXd x(df); 

	//initialize x with values in X. fill Qx with zeros. 
	for (int i = 0; i < df; i++) {
		x(i) = X[i];
	}
	Qx.fill(0.0);

	//fill Qx (Note: Qx is df by b after returning)
	makeQx(N, M, Qx, x);

	//get the orthogonal complement of Qx via the last d cols of Q (QR decomp of Qx)
	Eigen::MatrixXd Q(df,d); Q.fill(0.0);
	QRortho(Qx, Q, d); 

	//generate a proposal move, v, in the tangent space using Q
	Eigen::VectorXd v(df); v.fill(0.0);
	proposeTan(Q, v, d);

	//project the point to M, by adding w=Qx*a, to get sample, y
	Eigen::VectorXd a(b); a.fill(0.0);
	project(x+v, Qx);
}














}