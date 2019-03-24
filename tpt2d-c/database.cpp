#include "database.h"
#include <fstream>
#include <cstdlib>
namespace bd {


//state constructor
State::State() {
	am = NULL; coordinates = NULL; P = NULL;
	freq = 0; bond = 0; num = 0; denom = 0;
	num_coords = 0;
}

//state deconstructor
State::~State() {
	delete []am; delete []coordinates; delete []P;
}

//database constructor
Database::Database(int N_, int num_states_) {
	N = N_; num_states = num_states_;
	states = new State[num_states];
}

//database deconstructor
Database::~Database() {
	delete []states;
}

//sum the entries of s.P
int State::sumP(int num_states) const{
	int S = 0;
	for (int i = 0; i < num_states; i++) {
		S += P[i];
	}
	return S;
}

//pull a random set of coordinates from the available
const Cluster& State::getRandomIC() const {//todo fix rng
	int rand_state = rand() % num_coords;
	return coordinates[rand_state];
}


//function to read in the database and store in database class
Database* readData(std::string& filename) {
	std::ifstream in_str(filename);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		return NULL;
	}

	//read first line, N = number of particles
	int N;
	in_str >> N;

	//read second line - number of states
	int num_lines;
	in_str >> num_lines;


	bool val; //check if there is another line
	int coord; //number of sample coordinates for a state
	int index; //loop over states
	double x, y; //coordinates in a point
	char extra; //flag for whether mfpt estimates are in file

	//call the database class constructor
	Database* database = new Database(N, num_lines);

	//fill the database state classes
	while (in_str >> val) {
		//create reference to state, database[index]
		State& s = (*database)[index];

		//fill in adjacency matrix
		s.am = new bool[N*N];
		s.am[0] = val;
		for (int i = 1; i < N*N; i++) {
			in_str >> s.am[i];
		}

		//fill in frequency, bonds, and coords
		in_str >> s.freq;
		in_str >> s.bond;
		in_str >> coord;
		s.num_coords = coord;

		//fill in the sample coordinates
		s.coordinates = new Cluster[coord];
		for (int i = 0; i < coord; i++) {
			s.coordinates[i].setNumPoints(N);
			for (int j = 0; j < N; j++) {
				in_str >> x >> y;
				s.coordinates[i][j] = Point(x,y);
			}
		}

		//check the extra value for existence of mfpt estimates
		in_str >> extra;
		if (extra == 'N') {//no mfpt estimates, initialize to 0
			s.num = 0;
			s.denom = 0;
			s.P = new int[num_lines];
			for (int i = 0; i < num_lines; i ++) {
				s.P[i] = 0;
			}
		}
		else if (extra == 'Y') {//mfpt estimates exist, read in
			in_str >> s.num;
			in_str >> s.denom;
			s.P = new int[num_lines];
			for (int i = 0; i < num_lines; i ++) {
				in_str >> s.P[i];
			}
		}

		//next state
		index++;
	}
	in_str.close();
	return database;
}

//write function to output the updated database to a file
std::ostream& State::print(std::ostream& out_str, int N) const {
	for (int i = 0; i < N*N; i++) {
		out_str << am[i] << ' ';
	}
	out_str << freq << ' ';
	out_str << bond << ' ';
	out_str << num_coords << ' ';
	for (int i = 0; i < num_coords; i++) {
		for (int j = 0; j < N; j++) {
			out_str << coordinates[i][j].x << ' ' << coordinates[i][j].y << ' ';
		}
	}
	out_str << "N\n";
}

std::ostream& operator<<(std::ostream& out_str, const Database& db) {
	out_str << db.N << '\n';
	out_str << db.num_states << '\n';
	for (int i = 0; i < db.num_states; i++) {
		db[i].print(out_str, db.N);
	}
}

}