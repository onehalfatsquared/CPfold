#pragma once
#include <vector>
#include <ios>
namespace bd { 
class Database; 

class HCC;

/* HydroCluster class will hold the time series of a clusters positions, 
	 as well as some helpful auxiliary variables */

class HydroCluster {
	public:
		HydroCluster(); 
		~HydroCluster();
		

		double* clusters;
		int cNum;
		int N; int maxT;

	private:	

		//copy constructors - restricts compiling when user tries to copy a database
		HydroCluster(const HydroCluster&) {
			throw 1;
		}
		HydroCluster& operator=(const HydroCluster&) {
			throw 1;
		}
};

/* The HCC will be a collection of HydroClusters */

class HCC {
	public:
		HCC(int num_clusters_, int N_, int maxT_); 
		~HCC();
		HydroCluster& operator[](int index) {return hc[index];}
		const HydroCluster& operator[](int index) const {return hc[index];}
		//friend std::ostream& operator<<(std::ostream&, const Database&);

		HydroCluster* hc; 
		int N; int maxT; 

		//accessor functions
		int getN() const {return N;}
		int getMaxT() const {return maxT;}
		int getNumClusters() const {return num_clusters;}

	private:	
		int num_clusters; 

		//copy constructors - restricts compiling when user tries to copy a database
		HCC(const HCC&) {
			throw 1;
		}
		HCC& operator=(const HCC&) {
			throw 1;
		}
};

//input read functions
HCC* extractData(std::string& filename, int N, int maxT);

//transition detection functions
void getCoordinates(HydroCluster& hc, double* X, int N, int time);
void getAdj(double* X, int N, int* M, double cutoff);
void checkState(int N, double* X, int state, Database* db,
							  bool& reset, int& new_state);
void determineTransitions(HCC* hc, Database* db, double tps);
void determineTransitionTimes(HCC* hc, Database* db, double tps, std::ostream& ofile);
void determineTransitionStates(HCC* hc, Database* db, double tps, std::ostream& ofile);

//functions to compute statistics
void distributionFHT(HCC* hc, Database* db, std::vector<double>& q, int which);
void distributionFHT2(HCC* hc, Database* db, std::vector<double>& q, int which);
void timeAverageFHT(HCC* hc, Database* db, std::vector<double>& q, int which);
void sampleStats(std::vector<double> X, double& M, double& V);
void clustersFHT(HCC* hc, Database* db, std::ostream& ofile);

//test functions
void testExtract(HCC* hc);


}