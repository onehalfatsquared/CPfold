#include <math.h>
#include <stdio.h>
#include <string.h>
#include "database.h"
#include "bDynamics.h"
#include "tpt.h"
namespace bd{




void computeCommittor(double* q, double* T, int num_states, int initial, int* target) {
	//set up and solve equation for the committor function





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
			//fill in T_ji //todo
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
		mfpt = (*db)[state].getMFPT();
		S = (*db)[state].sumP(num_states); //get normalizing constant for this row

		//loop over row. Copy data into rate matrix
		for (int i = 0; i < num_states; i++) {
			T[toIndex(state, i, num_states)] = ((*db)[state].getP(i) / S) / mfpt;
		}
	}
}


void performTPT(int N, int initial, int target, Database* db) {
	//perform tpt calculations from initial to target states

	//construct the transition rate matrix 
	//step 1 - initialize the matrix
	int num_states = db->getNumStates();
	double* T = new double[num_states*num_states];
	for (int i = 0; i < num_states*num_states; i++) T[i] = 0;

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

	int* targets; //if there are multiple target states

	//solve dirichlet problem for q
	computeCommittor(q, T, num_states, initial, targets);



	//free the memory
	delete []T; delete []targets;

}




}