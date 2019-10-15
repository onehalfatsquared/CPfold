#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "point.h"
#include "database.h"
#include "bDynamics.h"
#include "tpt.h"
#include "nauty.h"
#include "graph.h"
#include "graphviz.h"
#include "../defines.h"
namespace bd{

void computeFlux(int num_states, double* q, double* T, double* eq, double* flux) {
	//compute the probability flux from generator, invariant measure, and committor

	int i, j;

	for (int entry = 0; entry < num_states*num_states; entry++) {
		index2ij(entry, num_states, i, j);
		if (i != j) {
			flux[entry] = eq[i]*T[entry]*q[j]*(1-q[i]);
		}
	}
}

void computeCommittor(double* q, double* T, int num_states, int initial, std::vector<int> targets) {
	//set up and solve equation for the committor function

	//initialize the matrix and vectors
	Eigen::MatrixXd TM(num_states,num_states); 
	Eigen::VectorXd b(num_states); 
	Eigen::VectorXd x(num_states);      


	//fill the matrix and RHS appropriately
	//copy T into TM, fill b with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		TM(i) = T[i]; 
	}
	b.fill(0.0);

	//loop over initial and target states. edit corresponding rows.
	//initial
	for (int j = 0; j < num_states; j++) {
		if (j != initial) {
			TM(initial, j) = 0;
		}
		else {
			TM(initial, j) = 1;
		}
	}

	//target
	for (int ind = 0; ind < targets.size(); ind++) {
		int target = targets[ind];
		for (int j = 0; j < num_states; j++) {
			if (j != target) {
				TM(target, j) = 0;
			}
			else {
				TM(target, j) = 1;
			}
		}
		b(target) = 1;
	}

	//solve the linear system
	x = TM.lu().solve(b);

	//store solution in q
	for (int i = 0; i < x.size(); i++) {
		q[i] = abs(x(i));
	}
}

void computeHittingProbability(double* P, int num_states, std::vector<int> endStates, 
															 double* U) {
	//compute the hitting probabilities for the states in endStates

	//init and fill all the required matrices
	Eigen::MatrixXd R(num_states,num_states); Eigen::MatrixXd Q(num_states,num_states); 
	Eigen::MatrixXd Pm(num_states,num_states); Eigen::MatrixXd Um(num_states,num_states); 
	for (int i = 0; i < num_states*num_states; i++) {
		R(i) = 0; Um(i) = 0; 
		Q(i) = P[i]; Pm(i) = P[i];
	}


	//let R be matrix with only columns corresponding to end states, diagonal = 1
	//let Q be matrix with rows/cols corresponding to end state zero'd out
	for (int i = 0; i < endStates.size(); i++) {
		int state = endStates[i];
		R.col(state) = Pm.col(state); R(state,state) = 1;
		Q.row(state) = Um.row(state); Q.col(state) = Um.col(state); 
	}

	//solve for Um via (I-Q)Um = R
	Eigen::MatrixXd D = Eigen::MatrixXd::Identity(num_states,num_states) - Q;
	Um = D.lu().solve(R);

	//store back in U array
	for (int i = 0; i < num_states*num_states; i++) {
		U[i] = Um(i);
	}
}


void fillDiag(double* T, int num_states) {
	//fill in diagonal elements with negative sum of entries

	double S;
	for (int state = 0; state < num_states; state++) {
		S = 0;
		for (int target = 0; target < num_states; target++) {
			S += T[toIndex(state, target, num_states)];
		}
		T[toIndex(state, state, num_states)] = -S;
	}
}


void satisfyDB(double* T, int num_states, Database* db, double* eq) {
	//make T satisfy detailed balance

	int row = 0; int column = 0; //row and column indices in matrix

	//if an entry is non-zero, use DB to fill its transpose
	for (int entry = 0; entry < num_states*num_states; entry++) {
		index2ij(entry, num_states, row, column);
		if (T[entry] > 0) {
			T[toIndex(column, row, num_states)] = T[entry] * eq[row] / eq[column];
		}
	}
}

void computeFreeEnergy(int num_states, double* Z, double* F) {
	//compute the free energy from configurational partition function

	double minF = 1e10;

	//compute the helmholtz free energy of each cluster from Z
	for (int i = 0; i < num_states; i++) {
		F[i] = -1.0/BETA * log(Z[i]);
		if (F[i] < minF) {
			minF = F[i];
		}
	}

	//shift to make 0 the free energy minimum
	for (int i = 0; i < num_states; i++) {
		F[i] -= minF;
	}
}
void computePartitionFn(int num_states, Database* db, double* Z) {
	//compute the configurational part of the partition fn from each cluster

	createMeasure(num_states, db, Z, KAP);
}

void createMeasure(int num_states, Database* db, double* eq, double kappa) {
	//create the equilibrium distribution for this problem

	double kap0 = KAP;                 //sticky parameter for initial measurement
	double stickyRatio = kappa / kap0; //ratio of current to old sticky parameters
	double beta = BETA;                //inverse temp

	double Z = 0;                     //normalizing constant for eq

	//loop over states, get the equilibrium measure entry
	for (int i = 0; i < num_states; i++) {
		//get initial eq prob and number of bonds
		int b = (*db)[i].getBonds();
		double prob = (*db)[i].getFrequency();

		//new eq prob is prob*(ratio)^(beta*bonds)
		double new_prob = prob*pow(stickyRatio,beta*b);

		//increment Z and add to array
		Z += new_prob;
		eq[i] = new_prob;
	}

	//re-normalize
	for (int i = 0; i < num_states; i++) eq[i] /= Z;
}

void createProbabilityMatrix(double* T, int num_states, double* P) {
	//create probability matrix from rate matrix

	int i,j; 

	for (int entry = 0; entry < num_states*num_states; entry++) {
		index2ij(entry, num_states, i, j);
		if (i != j) {
			P[entry] = - T[entry] / T[toIndex(i,i,num_states)];
		}
	}
}

void createTransitionMatrix(double* T, int num_states, Database* db, 
														std::vector<int>& endStates) {
	//create rate matrix from data in DB - forward rates

	//declare storage for each state
	double mfpt; int S;

	//loop over all states information
	for (int state = 0; state < num_states; state++) {
		//get all relevant data
		mfpt = (*db)[state].getMFPT();
		S = (*db)[state].sumP(); //get normalizing constant for this row

		//if S = 0, this is end state, add it to vector
		if (S == 0) {
			endStates.push_back(state);
		}

		//fill in value in transition matrix
		std::vector<Pair> P = (*db)[state].getP();
		for (int i = 0; i < P.size(); i++) {
			T[toIndex(state, P[i].index, num_states)] = (P[i].value / S) / mfpt;
		}
	}
}


void performTPT(int N, int initial, int target, Database* db, bool getIso) {
	//perform tpt calculations from initial to target states

	//parameters
	double kappa = 50.0;
	int num_states = db->getNumStates();
	std::vector<int> endStates;

	//init and compute configurational partition function and free energy
	double* Z = new double[num_states]; double* F = new double[num_states];
	for (int i = 0; i < num_states; i++) Z[i] = F[i] = 0;
	computePartitionFn(num_states, db, Z); 
	computeFreeEnergy(num_states, Z, F);

	//step 1 - construct rate matrix and probability transition matrix
	double* T = new double[num_states*num_states];
	double* P = new double[num_states*num_states];
	double* U = new double[num_states*num_states]; //hitting probability matrix
	for (int i = 0; i < num_states*num_states; i++) {
		T[i] = P[i] = U[i] = 0;
	}

	//step 2 - get bonds->bonds+1 entries from mfpt estimates
	createTransitionMatrix(T, num_states, db, endStates);

	//step 2.5 - make array for equilibrium distribution
	double* eq = new double[num_states];
	createMeasure(num_states, db, eq, kappa);

	//step 3 - fill in transposed entries such that T satisfies detailed balance
	satisfyDB(T, num_states, db, eq);

	//step 4 - fill in diagonal with negative sum of entries
	fillDiag(T, num_states);

	//step 5 - use the filled rate matrix to compute probability matrix
	createProbabilityMatrix(T, num_states, P);
	
	//solve for hitting probabilities to endStates states
	computeHittingProbability(P, num_states, endStates, U);

	//solve for the committor
	//initialize committor in q
	double* q = new double[num_states];
	for (int i = 0; i < num_states; i++) q[i] = 0;

	//build target state vector. check if including isomorphic states
	std::vector<int> targets; 
	if (getIso == 1) {
		//find all states isomorphic to target
		findIsomorphic(N, num_states, target, db, targets);
	}
	else if (getIso == 0) {
		targets.push_back(target);
	}

	/*
	//print the list of target states to check if correct
	for (int i = 0; i < targets.size(); i++) {
		printf("%d\n", targets[i]);
	}
	*/
	
	//solve dirichlet problem for committor, q 
	computeCommittor(q, T, num_states, initial, targets);

	//for (int i = 0; i < num_states; i++) std::cout << q[i] << "\n";


	//init and compute the probability fluxes
	double* flux = new double[num_states*num_states]; 
	for (int i = 0; i < num_states*num_states; i++) flux[i] = 0;
	computeFlux(num_states, q, T, eq, flux);

	//for (int i = 0; i < num_states*num_states; i++) std::cout << flux[i] << "\n";

	//make a graph 
	Graph* g = makeGraph(db);

	//print out graphviz
	printGraphRev(g, initial, F, flux, 1, 1, 0);

	//free the memory
	delete []T; delete []q; delete []eq; delete []flux;
	delete []Z; delete []F; delete []P ; delete [] U;

}




}