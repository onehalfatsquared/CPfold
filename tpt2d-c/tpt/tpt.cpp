#include <math.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <../Eigen/Dense>
#include "point.h"
#include "database.h"
#include "bDynamics.h"
#include "tpt.h"
#include "nauty.h"
namespace bd{



void computeCommittor(double* q, double* T, int num_states, int initial, std::vector<int> targets) {
	//set up and solve equation for the committor function

	//ex
	//num_states = 4; initial = 0; targets[0]=3;

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

	//ex
	/*TM << -1,1,0,0,
				0.25,-1,0.75,0,
				0,0.25,-1,0.75,
				0,0,1,-1;
	*/

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

	//ex
	//for (int i = 0; i < x.size(); i++) printf("%f\n", x(i));
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


void satisfyDB(double* T, int num_states, Database* db) {
	//make T satisfy detailed balance

	int row = 0; int column = 0;

	for (int entry = 0; entry < num_states*num_states; entry++) {
		index2ij(entry, num_states, row, column);
		if (T[entry] > 0) {
			//fill in T_ji - needs partitions functions? //todo
		}
	}
}


void createTransitionMatrix(double* T, int num_states, Database* db) {
	//create rate matrix from data in DB

	//declare storage for each state
	int num; int den; 
	double mfpt; int S;

	//loop over all states information
	for (int state = 0; state < num_states; state++) {
		//get all relevant data
		num = (*db)[state].getNumerator(); 
		den = (*db)[state].getDenominator(); 
		//todo
		mfpt = (*db)[state].getMFPT(); mfpt = 1; //FOR TESTING ONLY. DELETE LATER. NO DATA YET.
		S = (*db)[state].sumP(num_states); //get normalizing constant for this row
		S = 1; //DELETE LATER.

		//loop over row. Copy data into rate matrix
		for (int i = 0; i < num_states; i++) {
			T[toIndex(state, i, num_states)] = ((*db)[state].getP(i) / S) / mfpt;
		}
	}
}


void performTPT(int N, int initial, int target, Database* db, bool getIso) {
	//perform tpt calculations from initial to target states

	//construct the transition rate matrix 
	//step 1 - initialize the matrix
	int num_states = db->getNumStates();
	double* T = new double[num_states*num_states];
	for (int i = 0; i < num_states*num_states; i++) {
		T[i] = 0;
	}

	//step 2 - get bonds->bonds+1 entries from mfpt estimates
	createTransitionMatrix(T, num_states, db);

	//step 3 - fill in transposed entries such that T satisfies detailed balance
	satisfyDB(T, num_states, db);

	//step 4 - fill in diagonal with negative sum of entries
	fillDiag(T, num_states);

	//solve for the committor
	//initialize committor in q
	double* q = new double[num_states];
	for (int i = 0; i < num_states; i++) q[i] = 0;

	//build target state vector. check if inclusing isomorphic states
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

	std::vector<int> v;
	getPurgeStates(db, v);

	printf("size of v: %d\n", v.size());

	for (size_t i = 0; i < v.size(); i++) printf("%d\n", v[i]);

	printf("test\n");

	//free the memory
	delete []T; delete []q;
	printf("test2\n");

}




}