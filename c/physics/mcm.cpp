#include <math.h>
#include <random>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <chrono>
#include "bDynamics.h"
#include "database.h"
#include "adjacency.h"
#include "sampling.h"
#include "../defines.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/QR>
#include <omp.h>
namespace mcm {

double getResidual(int N, int* M, int b, Eigen::VectorXd p) {
	//get the residual of the constraint eqs in 2 norm at point p

	Eigen::VectorXd q(b); q.fill(0.0);
	evalConstraint(N, M, q, p); 
	return q.norm();
}

bool project(int N, int* M, int b, Eigen::VectorXd z, Eigen::MatrixXd Qz, Eigen::VectorXd& a) {
	//project the point onto the manifold using newtons method

	a.fill(0.0); //initial guess of 0
	int itc = 0; //init an iteration counter
	double residual = 1e6;  //big initial residual

	//set up the newton solve
	Eigen::MatrixXd J(b,b); J.fill(0.0);
	Eigen::MatrixXd Q2(b,z.size()); Q2.fill(0.0);
	Eigen::VectorXd F(b); F.fill(0.0);
	Eigen::VectorXd da(a.size()); da.fill(0.0);

	while (residual > N_TOL) {

		//construct the system
		J.fill(0.0); Q2.fill(0.0); F.fill(0.0);
		makeJ(N, M, Q2, z+Qz*a);
		J = Q2*Qz;
		evalConstraint(N, M, F, z+Qz*a);

		//solve the system
		da = J.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-F);

		//update
		a += da; itc += 1;
		residual = getResidual(N,M,b,z+Qz*a);

		//debug line
		//std::cout << residual << "\n";

		//check if max iteration number has passed, return false is fo
		if (itc > N_ITER) {
			return false;
		}
	}

	return true;

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

void makeJ(int N, int* M, Eigen::MatrixXd& J, Eigen::VectorXd x) {
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

				J(count, DIMENSION*i) = 2*XD; J(count, DIMENSION*j) = -2*XD;
				J(count, DIMENSION*i+1) = 2*YD; J(count, DIMENSION*j+1) = -2*YD;
#if (DIMENSION == 3) 
				J(count, DIMENSION*i+2) = 2*ZD; J(count, DIMENSION*j+2) = -2*ZD;
#endif

				count +=1;
			}
		}
	}

}

bool checkInequality(int N, int* M, Eigen::VectorXd y) {
	//check if the sample y violates any of the inequality constraints

	double XD, YD;
	double ZD = 0;

	//loop over and construct system
	for (int i = 0; i <N; i ++) {
		for (int j = i+1; j < N; j++) {
			if (M[bd::toIndex(i,j,N)] == 0) {
				//compute distances
				XD = y(DIMENSION*i) - y(DIMENSION*j);
				YD = y(DIMENSION*i+1) - y(DIMENSION*j+1);
#if (DIMENSION == 3) 
				ZD = y(DIMENSION*i+2) - y(DIMENSION*j+2);
#endif 

				//check if distance is less than 1
				if (XD*XD + YD*YD + ZD*ZD < 1) {
					return false;
				}
			}
		}
	}

	return true;
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

void getSample(int N, int* M, int df, int b, int d, Eigen::VectorXd& x) {
	//use the MCM algorithm to generate a new sample in x

	//init the gradient of the constraint matrix
	Eigen::MatrixXd Qx(b,df);
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
	Eigen::VectorXd y(df); y.fill(0.0);
	bool success = project(N, M, b, x+v, Qx, a);
	//if the projection fails, return with X_(n+1) = X_n
	if (!success) {
		return;
	}
	//if success, y is the proposed sample
	y = x+v+Qx*a;

	//check for inequality violations. if success = false, return.
	success = checkInequality(N, M, y);
	if (!success) {
		return;
	}

	//next, we check the reverse step
}

void runSampler(int N, double* X, int* M) {
	//run the sampling algo

	//get the number of bonds in the configuration, get DoF
	int df = DIMENSION*N;           //num degrees of freedom in ambient system
	int b = getNumBonds(N,M);      //num bonds
	int d = df-b;                  //effective dimension of system

	//setup the vector with eigen 
	Eigen::VectorXd x(df); 

	//initialize x with values in X. 
	for (int i = 0; i < df; i++) {
		x(i) = X[i];
	}

	//get a new sample
	getSample(N, M, df, b, d, x);

	
}














}