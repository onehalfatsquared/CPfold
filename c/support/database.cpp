#include "database.h"
#include "nauty.h"
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <fstream>
#include <cstdlib>
namespace bd {


//state constructor
State::State() {
	am = NULL; coordinates = NULL; 
	//P = NULL; Z = NULL; Zerr = NULL;
	freq = 0; bond = 0; num = 0; denom = 0;
	num_coords = 0; mfpt = 0; 
	sigma = 0; num_neighbors = 0;
}

//state deconstructor
State::~State() {
	delete []am; delete []coordinates; //delete []P;
	//delete []Z; delete []Zerr;
}

//database constructor
Database::Database(int N_, int num_states_) {
	N = N_; num_states = num_states_;
	states = new State[num_states];
	lumpMap = new int[num_states];
	for (int i = 0; i < num_states; i++) {
		lumpMap[i] = i;
	}
}

//database deconstructor
Database::~Database() {
	delete []states; delete []lumpMap;
}

//sum the entries of s.P
int State::sumP() const{
	int S = 0;
	for (int i = 0; i < num_neighbors; i++) {
		S += P[i].value;
	}
	return S;
}

//pull a random set of coordinates from the available
const Cluster& State::getRandomIC() const {
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
	int index = 0; //loop over states
	double x, y; //coordinates in a point
	char extra; //flag for whether mfpt estimates are in file
	int v1; double v2; //storing index and info in a pair
	double v3; double v4;

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
			s.mfpt = 0;
			s.sigma = 0;
			s.num_neighbors = 0;
		}
		else if (extra == 'Y') {//mfpt estimates exist, read in
			in_str >> s.num;
			in_str >> s.denom;
			in_str >> s.mfpt;
			in_str >> s.sigma; 
			in_str >> s.num_neighbors;

			for (int i = 0; i < s.num_neighbors; i ++) {
				in_str >> v1; in_str >> v2; in_str >> v3; in_str >> v4;
				s.P.push_back(Pair(v1,v2));
				s.Z.push_back(Pair(v1,v3));
				s.Zerr.push_back(Pair(v1,v4));
			}
		}

		//next state
		index++; 
	}
	in_str.close();
	return database;
}

//write functions to output the updated database to a file
std::ostream& State::print(std::ostream& out_str, int N, int* lumpMap) const {
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
	out_str << "Y" << ' ';
	out_str << num << ' ';
	out_str << denom << ' ';
	out_str << mfpt << ' ';
	out_str << sigma << ' ';
	out_str << num_neighbors << ' ';
	for (int i = 0; i < num_neighbors; i ++) {
		out_str << lumpMap[P[i].index] << ' ' << P[i].value << ' ';
		out_str << Z[i].value << ' ' << Zerr[i].value << ' ';
	}
	
	out_str << '\n';
	
}

std::ostream& operator<<(std::ostream& out_str, const Database& db) {
	//overload << to print out the database

	//check if any states are being purged
	int num_purge = db.toPurge.size();

	//print the non-purged states
	out_str << db.N << '\n';
	out_str << db.num_states - num_purge << '\n';
	for (int i = 0; i < db.num_states; i++) {
		if (std::find(db.toPurge.begin(), db.toPurge.end(), i) == db.toPurge.end()) {
			db[i].print(out_str, db.N, db.lumpMap);
		}
	}
}







void makeNM(int N, int state, Database* db, Eigen::VectorXd x, Eigen::MatrixXd& J, Eigen::VectorXd& F) {
	//make the matrix and vector to solve for Newtons method

	double XD, YD;
	int count = 0;

	//loop over and construct system
	for (int i = 0; i <N; i ++) {
		for (int j = i+1; j < N; j++) {
			if ((*db)[state].isInteracting(i,j,N)) {
				XD = x(2*i) - x(2*j);
				YD = x(2*i+1) - x(2*j+1);
				F(count) = XD*XD + YD*YD -1.0;
				J(count, 2*i) = 2*XD; J(count, 2*j) = -2*XD;
				J(count, 2*i+1) = 2*YD; J(count, 2*j+1) = -2*YD;
				count +=1;
			}
		}
	}
}


bool checkPhysicalState(int N, int state, Database* db) {
	//check if a state is physical - use as start point in Newtons method

	//check if there are any example coordinates
	int coord = (*db)[state].getNumCoords();
	if (coord == 0) {
		//printf("No coordinates\n");
		return false;
	}

	//set parameters to newtons method
	int max_iter = 50;
	double tol = 1e-8;

	//get the number of bonds
	int b = (*db)[state].getBonds(); 

	//initialize the matrix and vectors
	Eigen::MatrixXd J(b,2*N); 
	Eigen::VectorXd F(b); 
	Eigen::VectorXd dx(2*N);  
	Eigen::VectorXd x(2*N); 

	//initialize x. fill others with zeros.
	const Cluster& xc = (*db)[state].getRandomIC();
	for (int i = 0; i < N; i++) {
		x(2*i) = xc.points[i].x; x(2*i+1) = xc.points[i].y;
	}
	dx.fill(0.0); F.fill(0.0); J.fill(0.0);

	//output starting cluster - debug
	//for (int i =0; i < 2*N; i++) printf("%f\n", x(i));
	//printf("\n");


	//fill F and J. do solve with svd decomp.
	int iter;
	for (iter = 0; iter < max_iter; iter++) {
		makeNM(N, state, db, x, J, F);
		dx = J.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-F);
		x = x+dx;
		if (dx.norm() < tol) {
			break;
		}
	}

	//print message if no convergence
	if (iter == max_iter) {
		printf("State %d, Newton's Method terminated with residual %f\n", state, dx.norm());
	}

	//output final cluster - debug
	//for (int i =0; i < 2*N; i++) printf("%f\n", x(i));

	//compare initial and post newton adjacency matrix

	//set keywords for nauty
	int M = SETWORDSNEEDED(N);

	//initialize the nauty graphs and parameters
	graph g1[N*M]; graph g2[N*M];   //the starting graphs
	graph cg1[N*M]; graph cg2[N*M]; //graphs with canonical labeling
	int lab1[N], lab2[N];           //label arrays
	int ptn[N], orbits[N];          //needed to call functions
	static DEFAULTOPTIONS_GRAPH(options); //defualt options
	statsblk stats;                 //nauty statistics of graph
	options.getcanon = TRUE;        //get canonical labeling
	options.defaultptn = TRUE;      //ignore any coloring of graph

	//initialize arrays
	for (int i = 0; i < N; i++) {
		lab1[i] = lab2[i] = ptn[i] = orbits[i] = 0;
	}

	//build nauty graph of the old state, get canonical labeling
	buildNautyGraph(N, M, state, db, g1);
	EMPTYGRAPH(cg1, M, N);
	densenauty(g1, lab1, ptn, orbits, &options, &stats, M, N, cg1);

	//get adjacnecy matrix of post newton state
	int* AM = new int[N*N]; for (int i = 0; i < N*N; i++) AM[i] = 0;
	double* X = new double[2*N]; for (int i = 0; i < 2*N; i++) X[i]=x(i);
	getAdj(X, N, AM);

	//make graph
	//zero out any pre-existing graph
	EMPTYGRAPH(g2, M, N); 

	//add edges for every non-zero entry in adjacnecy matrix
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			if (AM[toIndex(i,j,N)] == 1) {
				ADDONEEDGE(g2, i, j, M);
			}
		}
	} 

	//free memory
	delete []X; delete []AM;

	//build nauty graph of the old state, get canonical labeling
	EMPTYGRAPH(cg2, M, N);
	densenauty(g2, lab2, ptn, orbits, &options, &stats, M, N, cg2);

	return checkIsomorphic(N, M, cg1, cg2);
}


void getPurgeStates(Database* db, std::vector<int>& toPurge) {
	//loop over all states. if unphysical, add to a list to purge.

	int N = db->getN(); int ns = db->getNumStates();

	for (int i = 0; i < ns; i++) {
		if (!checkPhysicalState(N, i, db)) {
			//state is unphysical, but check again to ensure IC not bad
			if (!checkPhysicalState(N, i, db)) {
				//state is unphysical, but check again to ensure IC not bad
				if (!checkPhysicalState(N, i, db)) {
					//state is unphysical. tell user and add to vector
					printf("State %d is unphysical.\n", i);
					toPurge.push_back(i);
				}
			}
		}
	}
}

void purgeUnphysical(Database* db) {
	//write to new file, skip purge states

	//init purge vector
	std::vector<int> toPurge;

	//get states to purge
	getPurgeStates(db, toPurge);

	//give to database
	db->toPurge = toPurge;

}

void combinePairs(std::vector<Pair>& p1, std::vector<Pair> p2) {
	//combines two vectors of pairs into one vector of pairs, p1.
	//if index is duplicated, the values are summed. 

	int i,j;

	bool flag;

	for ( i = 0; i < p2.size(); i++) {
		int index = p2[i].index; 
		flag = false;
		for ( j = 0; j < p1.size(); j++) {
			if (index == p1[j].index) {//both vectors have this state
				p1[j].value += p2[i].value;
				flag = true;
			}
		}
		if (!flag) {//only vector 2 has this, add to vector 1
			p1.push_back(p2[i]);
		}
	}
}

void getMinIndex(int N, int ns, std::vector<Pair> P, Database* db, 
									std::vector<int>& iso, std::vector<Pair>& isoPair ) {

	isoPair.clear();

	for (int j = 0; j < P.size(); j++) {
			iso.clear();
			findIsomorphic(N, ns, P[j].index, db, iso);
			int min_iso = *std::min_element(iso.begin(), iso.end());
			isoPair.push_back(Pair(min_iso, P[j].value));
			//todo - include Z and Zerr
		}
}

void lumpEntries(Database* db, int state, std::vector<int> perms) {
	//updates the entries of the lumped state //todo Z and Zerr

	int N = db->getN(); int ns = db->getNumStates();

	//basic quantities
	int freq = (*db)[state].getFrequency();
	int num = (*db)[state].getNumerator();
	int denom = (*db)[state].getDenominator();
	double mfpt = (*db)[state].getMFPT();
	double sigma = (*db)[state].getSigma();
	double dt = mfpt * denom / num;

	//state dependent quantities
	std::vector<Pair> P = (*db)[state].getP();
	std::vector<Pair> Z = (*db)[state].getZ();
	std::vector<Pair> Zerr = (*db)[state].getZerr();

	//isomorphism vectors
	std::vector<int> iso;
	std::vector<Pair> isoPair;

	//make new P that contains min isomorphism indexes
	std::vector<Pair> Pnew; Pnew.clear();

	for (int i = 0; i < perms.size(); i++) {

		printf("%d, ", perms[i]);
		//update basic quantities
		freq += (*db)[perms[i]].getFrequency();
		num += (*db)[perms[i]].getNumerator();
		denom += (*db)[perms[i]].getDenominator();
		double nextMFPT = (*db)[perms[i]].getMFPT();
		double nextSQ = nextMFPT*nextMFPT;
		double nextSigma = (*db)[perms[i]].getSigma();
		//progogate error in mfpt reciprocals
		mfpt = 1 / (1/mfpt + 1/nextMFPT); 
		sigma = sqrt(sigma*sigma + nextSigma*nextSigma/(nextSQ*nextSQ));

		//update state dependent quantities
		//get quantities for states getting lumped
		std::vector<Pair> P2 = (*db)[perms[i]].getP();
		std::vector<Pair> Z2 = (*db)[perms[i]].getZ();
		std::vector<Pair> Zerr2 = (*db)[perms[i]].getZerr();

		//compute the minimum index of this states isomorphisms
		getMinIndex(N, ns, P2, db, iso, isoPair);

		//have the pair vector to update P with. call function to update
		combinePairs(Pnew, isoPair);


	}

	//write the updates back to the DB
	(*db)[state].freq = freq;
	(*db)[state].num = num;
	(*db)[state].denom = denom;
	(*db)[state].mfpt = mfpt;
	(*db)[state].sigma = sigma*mfpt*mfpt;
	(*db)[state].num_neighbors = Pnew.size();
	(*db)[state].P = Pnew;
	(*db)[state].Z = Z; // todo
	(*db)[state].Zerr = Zerr;// todo
		

}


void lumpPerms(Database* db) {
	//loop over states, lump perms, get vector of repeated states

	int N = db->getN(); int ns = db->getNumStates();
	std::vector<int> repeated;
	std::vector<int> perms;
	int* lumpMap = new int[ns]; 

	int count = 0;

	for (int i = 0; i < ns; i++) {
		//check that i is not in repeated list
		if (std::find(repeated.begin(), repeated.end(), i) == repeated.end()) {
			//find iso states, lump them together
			findIsomorphic(N, ns, i, db, perms);
			printf("State %d is a permutation of state ", i);
			lumpEntries(db, i, perms);
			printf("\n");
			//add to repeated vector, so they are not called twice
			repeated.insert(repeated.end(), perms.begin()+1, perms.end());
			//update lumpMap
			for (int j = 0; j < perms.size(); j++) {
				lumpMap[perms[j]] = count;
			}
			count += 1;
			//clear perms vector
			perms.clear();
		}
	}
	//store the purge vector
	db->toPurge = repeated;
	
	//store lumpMap in DB
	for (int i = 0; i < ns; i++) {
		db->lumpMap[i] = lumpMap[i];
	}

	//free memory
	delete []lumpMap;

}

void combineMFPTdata(Database* db1, Database* db2) {
	//combine the mfpt data in the two databases, store in the first
	
	//get database properties, assumes both have same num states
	int N = db1->getN(); int ns = db1->getNumStates();

	//loop over the states and combine estimates
	for (int state = 0; state < ns; state++) {
		//get new mfpt and error bar estimate - minimize var of a linear combination
		double mfpt1 = (*db1)[state].getMFPT(); double mfpt2 = (*db2)[state].getMFPT();
		double sigma1 = (*db1)[state].getSigma(); double sigma2 = (*db2)[state].getSigma();
		double v1 = sigma1*sigma1; double v2 = sigma2*sigma2;
		double S = 1.0/v1 + 1.0/v2;
		double mfpt = 1.0/S * (mfpt1/v1 + mfpt2/v2);
		double sigma = 1.0/S * sqrt(1.0/v1 + 1.0/v2);

		//combine the hit counts via pairs
		std::vector<Pair> p1 = (*db1)[state].getP();
		std::vector<Pair> p2 = (*db2)[state].getP();
		combinePairs(p1, p2);

		//add back to database 1
		(*db1)[state].mfpt = mfpt;
		(*db1)[state].sigma = sigma;
		(*db1)[state].num_neighbors = p1.size();
		(*db1)[state].P = p1;
	}
}

void printRatios(Database* db1, Database* db2) {
	//print the ratio of mfpt estimates in two databases to compare

	//get database properties, assumes both have same num states
	int N = db1->getN(); int ns = db1->getNumStates();

	//loop over the states and combine estimates
	for (int state = 0; state < ns; state++) {
		//get mfpt estimates
		double mfpt1 = (*db1)[state].getMFPT(); 
		double mfpt2 = (*db2)[state].getMFPT();
		double ratio = mfpt1/mfpt2;

		printf("State %d has an mfpt ratio %f\n", state, ratio);
	}
}

void printProbs(Database* db1, Database* db2) {
	//print the probability distributions for each state in two dbs

	//get database properties, assumes both have same num states
	int N = db1->getN(); int ns = db1->getNumStates();

	//loop over the states and combine estimates
	for (int state = 0; state < ns; state++) {
		//get mfpt estimates
		int S1 = (*db1)[state].sumP();
		int S2 = (*db2)[state].sumP();

		std::vector<Pair> p1 = (*db1)[state].getP();
		std::vector<Pair> p2 = (*db2)[state].getP();

		printf("State %d transition probabilities\n", state);

		for (int i = 0; i < p1.size(); i++) {
			int index = p1[i].index; int val1 = p1[i].value;
			int val2 = -1;
			for (int j = 0; j < p2.size(); j++) {
				if (p2[j].index == index) {
					val2 = p2[j].value;
					break;
				}
			}
			double prob1 = float(val1)/float(S1);
			double prob2 = float(val2)/float(S2);
			double diff = abs(prob2-prob1);

			printf("%d, %f, %f, Difference = %f\n", index, prob1, prob2, diff);
		}
	}

}

void updateFreq(Database* db, std::string filename) {
	//update frequency value with eq prob from better simulation
	
}



}