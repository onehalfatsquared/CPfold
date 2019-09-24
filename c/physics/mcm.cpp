#include <math.h>
#include <random>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
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

/****************************************************/
/*********** Code for Trimer Test *******************/
/****************************************************/


double getBondAngle(Eigen::VectorXd x) {
	//get the bond angle for a trimer

	//init the storage for vectors between particles (1,2) and (2,3)
	Eigen::VectorXd v1(DIMENSION); v1.fill(0.0);
	Eigen::VectorXd v2(DIMENSION); v2.fill(0.0);

	//fill the displacement vectors
	for (int i = 0; i < DIMENSION; i++) {
		v1(i) = x(DIMENSION+i)-x(i);
		v2(i) = x(2*DIMENSION+i)-x(DIMENSION+i);
	}

	//compute the scaled dot product
	double sdp = (v1.dot(v2))/(v1.norm()*v2.norm());

	//return the inverse cosine of sdp to get angle
	return acos(sdp);
}

void runTestTrimer(int N, double* X, int* M, int num_samples) {
	//run the sampling algo on a trimer test problem

	//get the number of bonds in the configuration, get DoF
	int df = DIMENSION*N;           //num degrees of freedom in ambient system
	int b = getNumBonds(N,M);      //num bonds
	int d = df-b;                  //effective dimension of system

	//setup the vector with eigen 
	Eigen::VectorXd x(df);        //store the current sample
	Eigen::VectorXd prev(df);     //stores the previous sample 

	//initialize x with values in X. 
	for (int i = 0; i < df; i++) {
		x(i) = X[i];
	}

	//construct an rng class
	RandomNo* rngee = new RandomNo();

	//init a counter
	int ctr = 0; int burnIn = 1000;

	//output file
	std::ofstream ofile;
	ofile.open("ba.txt");

	//get some samples
	while (ctr < num_samples+burnIn) {
		prev = x;
		getSample(N, M, df, b, d, x, rngee);
		if (ctr>burnIn)
			ofile << getBondAngle(x) << "\n";
		ctr++;
	}

	//close the file
	ofile.close();

	//free the memory
	delete rngee;
}

/****************************************************/
/*********** End code for trimer test ***************/
/****************************************************/


/****************************************************/
/*********** MCMC Internal Functions ****************/
/****************************************************/


double getResidual(int N, int* M, int b, Eigen::VectorXd p) {
	//get the residual of the constraint eqs in 2 norm at point p

	Eigen::VectorXd q(b); q.fill(0.0);
	evalConstraint(N, M, q, p); 
	return q.norm();
}

double getMH(Eigen::MatrixXd Qx, Eigen::MatrixXd Qy, Eigen::VectorXd x, 
						 Eigen::VectorXd y, Eigen::VectorXd v, Eigen::VectorXd vr, int d, int b) {
	//compute the Metropolis-Hastings acceptance probability for the move

	double N = fEval(b, Qy)*evalDensity(vr,d);
	double D = fEval(b, Qx)*evalDensity(v,d);
	double p = N/D;

	if (p > 1) 
		return 1.0;
	else
		return p;
}

double getJacobian(Eigen::MatrixXd Q) {
	//evaluate the delta fn jacobian term, inverse of pseudo-determinant of Q

	Eigen::VectorXd s = Q.jacobiSvd().singularValues();
  double pdet = 1;
  for (int i = 0; i < s.size(); i++) {
    if (abs(s(i)) > N_TOL) {  
      pdet *= abs(s(i));
    }
  }

  return 1.0/pdet;
}

double fEval(int b, Eigen::MatrixXd Q) {
	//evaluate the function to be sampled 

	double f = 1.0;
	for (int i = 0; i < b; i++) {
		f *= KAP;
	}

	return f*getJacobian(Q);
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

	//perform iterations of NM
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

	/* //Normalization not needed for MCMC
	for (int i = 0; i < d; i++) {
		Z /= (sqrt(2*M_PI)*SIG);
	}
	*/

	return Z*exp(-V/(2*SIG*SIG));
}

void proposeTan(Eigen::MatrixXd Q, Eigen::VectorXd& v, int d, RandomNo* rngee) {
	//generate a proposal move tangent to M - isotropic gaussian

	//store the weights of each vector
	Eigen::VectorXd r(d);

	//compute weights, standard dev, SIG, times N(0,1)
	for (int i = 0; i < d; i++) {
		r(i) = SIG*rngee->getG(); 
	}

	//get proposal move
	v = Q*r;
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
	//construct the jacobian matrix for use in newtons method

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

bool getSample(int N, int* M, int df, int b, int d, Eigen::VectorXd& x, RandomNo* rngee) {
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
	proposeTan(Q, v, d, rngee);

	//project the point to M, by adding w=Qx*a, to get sample, y
	Eigen::VectorXd a(b); a.fill(0.0);
	Eigen::VectorXd y(df); y.fill(0.0);
	bool success = project(N, M, b, x+v, Qx, a);
	//if the projection fails, return with X_(n+1) = X_n
	if (!success) {
		return false;
	}
	//if success, y is the proposed sample
	y = x+v+Qx*a;

	//check for inequality violations. if success = false, return.
	success = checkInequality(N, M, y);
	if (!success) {
		return false;
	}

	//next, we check the reverse step.
	// Get the matrix Qy
	Eigen::MatrixXd Qy(b,df);
	Qy.fill(0.0);
	makeQx(N, M, Qy, y);

	//get the orthogonal complement of Qy via the last d cols of Q (QR decomp of Qx)
	Q.fill(0.0);
	QRortho(Qy, Q, d); 

	//project x-y onto the space spanned by cols of Q to get reverse step, vr
	Eigen::VectorXd vr(df); vr.fill(0.0);
	for (int i = 0; i < d; i++ ) {
		Eigen::VectorXd u = Q.col(i);
		vr += u.dot(x-y)*u;
	}
	
	//compute the acceptance prob 
	double prob = getMH(Qx,Qy,x,y,v,vr,d,b);

	//accept or reject
	double U = rngee->getU();
	if (U > prob) {
		return false;
	}

	//if we reach here, y is accepted. Just need to check reverse projection
	success = project(N, M, b, y+vr, Qy, a);
	if (!success) {
		return false;
	}

	//projection successful, check if result gives back x
	Eigen::VectorXd xDiff(df); xDiff.fill(0.0);
	Eigen::VectorXd xr(df); xr.fill(0.0);
	xr = y+vr+Qy*a; xDiff = x-xr;
	if (xDiff.norm() < 100*N_TOL) {
		x=y;
		return true;
	}

	return false;
}

/****************************************************/
/*********** End MCMC Internal Functions ************/
/****************************************************/

/****************************************************/
/*********** MFPT Estimator Functions ***************/
/****************************************************/

void sampleStats(std::vector<double> X, double& M, double& V) {
	//return the sample mean of the vector X

	//init the mean and variance at 0, get sample size
	M = 0; V = 0;
	int N = X.size();

	//compute the mean
	for (int i = 0; i < N; i++) M += X[i];
	M /= float(N);

	//compute the variance 
	for (int i = 0; i < N; i++) V += (X[i]-M) * (X[i]-M);
	V /= (N-1);
}

void minVarEstimate(int sampleSize, double* means, double* variances) {
	//combine the estimated means in a linear combination to minimize variance


}

void equilibrate() {

}

void getSampleMFPT() {
	
}



void estimateMFPT(int N, int state, bd::Database* db) {
	/*estimate mean first passage time starting in state and going to state with
	one additional bond. Uses parallel implementations of a single walker with
	long trajectory.*/

	//set parameters
	int samples = SAMPLES; //number of hits per walker for estimator
	int eq = EQ; //number of steps to equilibrate for
	int num_states = db->getNumStates(); //total number of states

	//quantities to update - new estimates
	std::vector<bd::Pair> PM; std::vector<bd::Pair> PMshare;
	double mfpt = 0;
	double sigma = 0;

	//output start message for this state
	printf("Beginning MFPT Estimator for state %d out of %d.\n", state, num_states);

	//store mfpt estimates on each thread to get standard deviation
	double* mfptSamples; double* mfptVar; int num_threads;

	//debug line
	for (int i = 0; i < PM.size(); i++) printf("0 thread has:\n %d, %f\n", PM[i].index, PM[i].value);

	//open parallel region
	#pragma omp parallel private(PM) shared(PMshare) 
	{
		//initialize final samples storage - only on one processor - then barrier
		num_threads = omp_get_num_threads();
		#pragma omp single
		{
			mfptSamples = new double[num_threads];
			mfptVar     = new double[num_threads];
			for (int i = 0; i < num_threads; i++) {
				mfptSamples[i] = 0; mfptVar[i] = 0;
			}
		}
		#pragma omp barrier

		//init the private mfpt sample storage
		std::vector<double> mfptVec;

		//get starting coordinates randomly from the database
		const bd::Cluster& c = (*db)[state].getRandomIC();

		//cluster structs to arrays
		double* X = new double[DIMENSION*N];
		if (DIMENSION == 2)
			c.makeArray2d(X, N);
		else if (DIMENSION == 3) 
			c.makeArray3d(X, N);

		for (int step = 0; step < samples; step++) {
			//equilibrate the trajectories
			equilibrate();

			//get a sample - has to update PM
			getSampleMFPT();

			//add sample to mfptVec

			//revert to original coordinates

		}

		//get sample means and variances
		double M; double V;
		sampleStats(mfptVec, M, V);
		mfptSamples[omp_get_thread_num()] = M;
		mfptVar[omp_get_thread_num()] = V;

		//do update on PM vectors - need barrier - one at a time
		if (omp_get_thread_num() == 0) {
			PMshare = PM;
		}
		#pragma omp barrier
		#pragma omp critical
		{
			if (omp_get_thread_num() != 0) {
				combinePairs(PMshare, PM);
			}
		}

		//free cluster memory
		delete []X;

		//end parallel region
	}

	//update estimates
	//combine the mfptSamples entries to get min variance estimator
		

	//make a Z vector with same num of elements as P
	std::vector<bd::Pair> Z; 
	for (int i = 0; i < PMshare.size(); i++) {
		Z.push_back(bd::Pair(PMshare[i].index, 0));
	}

	//update database
	(*db)[state].mfpt = mfpt;
	(*db)[state].num = 0;
	(*db)[state].denom = 0;
	(*db)[state].num_neighbors = PMshare.size();
	(*db)[state].P = PMshare;
	(*db)[state].Z = Z;
	(*db)[state].Zerr = Z;
	(*db)[state].sigma = sigma;

	//print out final estimates - debug
	/*
	printf("%Numerator = %d, Denominator = %d, \n", num, den);
	double sum = 0;
	for (int i = 0; i < PMshare.size(); i++) {
		printf("State = %d, visits = %f\n", PMshare[i].index, PMshare[i].value);
		sum +=PMshare[i].value;
	}
	printf("sum of hits = %f\n", sum);
	printf("Total Estimate = %f +- %f\n", mfpt, sigma);
	for (int i = 0; i < num_threads; i++) printf("MFPT estimate %d = %f\n", i, mfptSamples[i]);
	*/


	//free memory
	delete []mfptSamples; delete []mfptVar;
}






}