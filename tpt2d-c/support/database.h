#pragma once
#include "point.h"
#include <string>

/*database structure to store all states.
	  N - number of particles in system
		num_states - number of known states

	state structure to store all info about a state.
		am - adjacency matrix for the state - N by N
		freq - the frequency of the state during a bd simulation
		bond - the number of bonds the state has
		coordinates - an array of sample coordinates in the state - unknown by 2*N;
		num - the numerator in an mfpt estimator
		den - the denominator in an mfpt estimator
		P - un-normalized transition probabilites out of this state - 1 by num_clusters (known)
*/

namespace bd {

class Database;

class State{
	public:
		State();
		~State();
		friend Database* readData(std::string& filename);

		//accessor functions
		int getFrequency() const {return freq;}
		int getBonds() const {return bond;}
		const Cluster& getRandomIC() const;
		bool isInteracting(int i, int j, int N) const {return am[i*N+j];}
		int getNumerator() const {return num;}
		int getDenominator() const {return denom;}
		int getP(int i) const {return P[i];}
		int sumP(int num_states) const;

		//print function
		std::ostream& print(std::ostream&, int) const;

	private:
		int freq, bond;
		int num, denom;
		bool* am; int* P;
		int num_coords; 
		Cluster* coordinates; 

		//copy constructors - restricts compiling when user tries to copy state
		State(const State&) {
			throw 1;
		}
		State& operator=(const State&) {
			throw 1;
		}


};

class Database {
	public:
		Database(int N_, int num_states_); 
		~Database();
		State& operator[](int index) {return states[index];}
		const State& operator[](int index) const {return states[index];}
		friend std::ostream& operator<<(std::ostream&, const Database&);

		//accessor functions
		int getN() const {return N;}
		int getNumStates() const {return num_states;}

	private:	
		int N; int num_states; State* states; 

		//copy constructors - restricts compiling when user tries to copy a database
		Database(const Database&) {
			throw 1;
		}
		Database& operator=(const Database&) {
			throw 1;
		}

};


//read/write data functions
Database* readData(std::string& filename);
std::ostream& operator<<(std::ostream&, const Database&);

}
