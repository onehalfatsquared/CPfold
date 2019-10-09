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
#include "../defines.h"
namespace bd{



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
		q[i] = x(i);
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

		std::cout << prob << "\n";

		//new eq prob is prob*(ratio)^(beta*bonds)
		double new_prob = prob*pow(stickyRatio,beta*b);

		//increment Z and add to array
		Z += new_prob;
		eq[i] = new_prob;
	}

	//re-normalize
	for (int i = 0; i < num_states; i++) eq[i] /= Z;
}


void createTransitionMatrix(double* T, int num_states, Database* db) {
	//create rate matrix from data in DB - forward rates

	//declare storage for each state
	double mfpt; int S;

	//loop over all states information
	for (int state = 0; state < num_states; state++) {
		//get all relevant data
		mfpt = (*db)[state].getMFPT();
		S = (*db)[state].sumP(); //get normalizing constant for this row

		//loop over row. Copy data into rate matrix
		/*
		for (int i = 0; i < num_states; i++) {
			T[toIndex(state, i, num_states)] = ((*db)[state].getP(i) / S) / mfpt;
		}
		*/
		std::vector<Pair> P = (*db)[state].getP();
		for (int i = 0; i < P.size(); i++) {
			T[toIndex(state, P[i].index, num_states)] = (P[i].value / S) / mfpt;
		}
	}
}


void performTPT(int N, int initial, int target, Database* db, bool getIso) {
	//perform tpt calculations from initial to target states

	//construct the transition rate matrix 

	double kappa = 4;

	//step 1 - initialize the matrix
	int num_states = db->getNumStates();
	double* T = new double[num_states*num_states];
	for (int i = 0; i < num_states*num_states; i++) {
		T[i] = 0;
	}

	//step 2 - get bonds->bonds+1 entries from mfpt estimates
	createTransitionMatrix(T, num_states, db);
	//for (int i = 0 ; i < num_states*num_states; i++) std::cout << T[i] << "\n";

	//step 2.5 - make array for equilibrium distribution
	double* eq = new double[num_states];
	createMeasure(num_states, db, eq, kappa);
	//for (int i = 0; i < num_states; i++) std::cout << eq[i] << "\n";

	//step 3 - fill in transposed entries such that T satisfies detailed balance
	satisfyDB(T, num_states, db, eq);

	//step 4 - fill in diagonal with negative sum of entries
	fillDiag(T, num_states);

	//for (int i = 0 ; i < num_states*num_states; i++) std::cout << T[i] << "\n";

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
	
	//solve dirichlet problem for q 
	computeCommittor(q, T, num_states, initial, targets);

	//free the memory
	delete []T; delete []q; delete []eq;

}




}