#pragma once
#include "point.h"
#include "pair.h"
#include "adjacency.h"
#include "nauty.h"
#include <../Eigen/Dense>
#include <vector>
#include <string>

/*database structure to store all states.
	  N - number of particles in system
		num_states - number of known states
		toPurge - states to not include when writing database to a new file
		lumpMap - a mapping of all original states to a new index upon lumping

	state structure to store all info about a state.
		am - adjacency matrix for the state - N by N
		freq - the frequency of the state during a bd simulation
		bond - the number of bonds the state has
		coordinates - an array of sample coordinates in the state - unknown by 2*N;
		num - the numerator in an mfpt estimator
		denom - the denominator in an mfpt estimator
		mfpt - num/den*delta_t
		sigma - standard deviation of the mfpt estimate
		num_neighbors - number of states that this state talks to
		P - un-normalized tranisition probabilities. stored sparsely in pairs
		Z - ratio of partition functions. sparse pair storage
		Zerr - standard deviation of partition function estimate

		values appear in this order when reading and writing to file.

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
		int getNumCoords() const {return num_coords;}
		const Cluster& getRandomIC() const;
		bool isInteracting(int i, int j, int N) const {return am[i*N+j];}
		int getNumerator() const {return num;}
		int getDenominator() const {return denom;}
		double getMFPT() const {return mfpt;}
		double getSigma() const {return sigma;}
		int sumP() const;
		std::vector<Pair> getP() const {return P;}
		std::vector<Pair> getZ() const {return Z;}
		std::vector<Pair> getZerr() const {return Zerr;}

		//quantities to update
		int freq, num, denom;
		int num_neighbors;
		double mfpt, sigma;
		//Pair* P; Pair* Z; Pair* Zerr;
		std::vector<Pair> P; std::vector<Pair> Z; std::vector<Pair> Zerr;

		//print function
		std::ostream& print(std::ostream&, int, int*) const;

	private:
		int bond;
		bool* am; 
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

		std::vector<int> toPurge;
		int* lumpMap;

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

bool checkPhysicalState(int N, int state, Database* db);
void makeNM(int N, int state, int b, Database* db, Eigen::VectorXd , Eigen::MatrixXd& , 
																							     Eigen::VectorXd& , Eigen::VectorXd);

void getPurgeStates(Database* db, std::vector<int>& toPurge);
void purgeUnphysical(Database* db);
void combinePairs(std::vector<Pair>& p1, std::vector<Pair> p2);
void lumpEntries(Database* db, int state, std::vector<int> perms);
void lumpPerms(Database* db);
void getMinIndex(int N, int ns, std::vector<Pair> P, Database* db, 
									std::vector<int>& iso, std::vector<Pair>& isoPair );

}
