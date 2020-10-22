#pragma once
#include "point.h"
#include "bDynamics.h"
#include "tpt.h"
#include <vector>
#include <map>
#include <deque>
#include <random>
#include "dlib/optimization.h"
namespace bd { 

typedef dlib::matrix<double,0,1> column_vector;

class Database; 

//create a class to store everything constant about a problem
class OptInfo {
	public:
		//constructors
		OptInfo();
		~OptInfo();

		OptInfo(int N_, int num_states_, int numTypes, int numInteractions_, int* particleTypes_, 
						int initial_, int target_, double rho_, Database* db_);
		OptInfo(int N, int num_states, int initial, int target, 
						double E, double rho, Database* db);
		OptInfo(int N, int initial, int target, 
						double h1, double h2);
		OptInfo(int N_, int num_states_, int initial_, int target_, 
						const Eigen::MatrixXd& h_);
		OptInfo(int N, int num_states, int initial, int target, 
						double h1, double h2);

		//call the optimization problem functions
		void CreateRateMatrix(Eigen::MatrixXd& R, double beta);
		void CreateRateMatrix(Eigen::MatrixXd& R, Eigen::MatrixXd& D, double beta);
		void CreateRateMatrix(Eigen::MatrixXd& R, const Eigen::MatrixXd& r_form, 
												  const Eigen::VectorXd& kappaVals);
		void CreateRateMatrix(Eigen::MatrixXd& R, Eigen::MatrixXd& D, 
													const Eigen::MatrixXd& r_form, const Eigen::VectorXd& kappa,
													int ixn);

		int GetProblemType() const {return problem_type;}
		int GetNumInteractions() const {return numInteractions;}




	private:
		int problem_type;
		//variables common to all problem types
		int N; int num_states;
		int initial; int target;

		//variables for problem type 0
		int numInteractions; int numTypes;
		int* particleTypes;

		//variables for problem type 1
		double E; double rho; 
		Database* db;

		//variables for problem type 2
		double h1; double h2;

		//variables for problem type 3
		Eigen::MatrixXd h;

};

inline OptInfo::OptInfo() {
	problem_type = 0;
	N = 0; num_states = 0; initial = 0; target = 0;
	E = 0; rho = 0;
	h1 = 0; h2 = 0;
	numInteractions = 0; numTypes = 0;
	particleTypes = NULL;
	db = NULL;
}

inline OptInfo::~OptInfo() {
	if (db != NULL) {
		delete db;
	}
	if (particleTypes != NULL) {
		delete []particleTypes;
	}
}

/****************************************************************/
/******************** Temperature Protocols *********************/
/****************************************************************/

//matrix creation
//self assembly problem
void createRateMatrix(Eigen::MatrixXd& R, int num_states, Database* db, 
												double beta, double E, double rho);
void createRateMatrix(Eigen::MatrixXd& R, Eigen::MatrixXd& D, int num_states, Database* db, 
												double beta, double E, double rho);

//2 state problem
void createRateMatrix(Eigen::MatrixXd& R, double beta, double h1, double h2);
void createRateMatrix(Eigen::MatrixXd& R, Eigen::MatrixXd& D, double beta,
											double h1, double h2);

//model problem
void createRateMatrix(Eigen::MatrixXd& R, int num_states, double beta, 
											const Eigen::MatrixXd& h);
void createRateMatrix(Eigen::MatrixXd& R, Eigen::MatrixXd& D, int num_states, double beta, 
											const Eigen::MatrixXd& h);

//1 big row problem
void createRateMatrix(Eigen::MatrixXd& R, int num_states, double beta, double h1, double h2);
void createRateMatrix(Eigen::MatrixXd& R, Eigen::MatrixXd& D, int num_states, double beta,
											double h1, double h2);


//helper functions
void applyLowerBound(Eigen::VectorXd& vec, int length, double bound);

void performGD(int N, int max_iters, Eigen::VectorXd& betaP, int num_states, 
							 OptInfo* problem, const Eigen::VectorXd& t_disc, 
							 const Eigen::VectorXd& t_centers, int initial, int target, double step,
							 double lowerBound, double grad_tol, double input_tol, bool verbose);

void createTransitionOperator(int N, int num_states, OptInfo* problem, 
															const Eigen::VectorXd& betaP, const Eigen::VectorXd& t_disc, 
															double t0, double tf, int interval, 
															Eigen::MatrixXd& trans_op);

double computeFinalProb(int N, int num_states, OptInfo* problem,
												const Eigen::VectorXd& betaP, const Eigen::VectorXd& t_disc,  
												int initial, int target);



//run the gradient descent
void findProtocol(int np, Database* db, int initial, int target);
void findProtocol2state();
void getGaps2state();
void flowerTestConstantBeta(int np, Database* db);

void modelFlowerSystem();

//testing
void testProtocol(int np, Database* db, int initial, int target);
void testRowSystemProtocols();

//dlib
void modelFlowerSystem_dlib();
void modelRowSystem_dlib();
void findProtocolSA_dlib(int np, Database* db, int initial, int target);


/****************************************************************/
/******************** Vector Valued Protocols *******************/
/****************************************************************/

void applyBounds(int elements, Eigen::VectorXd& kappaVec, 
								 double lower, double upper);

//matrix version
double computeFinalProb(int N, int num_states, OptInfo* problem, const Eigen::MatrixXd& f_rate,
												const Eigen::MatrixXd& kappaMat, const Eigen::VectorXd& t_disc,  
												int initial, std::vector<int> targets);
//vector version
double computeFinalProb(int N, int num_states, OptInfo* problem, const Eigen::MatrixXd& f_rate,
												const Eigen::VectorXd& kappaVec, const Eigen::VectorXd& t_disc,  
												int initial, std::vector<int> targets);

void evalGradient(int N, const Eigen::MatrixXd& kappaMat, int num_states, 
							 OptInfo* problem, const Eigen::MatrixXd& f_rate, const Eigen::VectorXd& t_disc, 
							 const Eigen::VectorXd& t_centers, int initial, std::vector<int> targets,
							 Eigen::VectorXd& R);

void createTransitionOperator(int N, int num_states, OptInfo* problem, 
															const Eigen::MatrixXd& f_rate, const Eigen::MatrixXd& kappaMat,
															const Eigen::VectorXd& t_disc, double t0, double tf, int interval, 
															Eigen::MatrixXd& trans_op);

double dlib_objective_call(int N, int num_states, int numInteractions, OptInfo* problem, 
													 column_vector& kappaVec, const Eigen::MatrixXd& f_rate, 
													 const Eigen::VectorXd& t_disc, const Eigen::VectorXd& t_centers, 
													 int initial, std::vector<int> targets);

void findProtocolSA_dlib_E(int np, Database* db, int initial, int target, bool useFile);
void findProtocolSA_adjoint(int np, Database* db, int initial, int target, bool useFile);
void testProtocol_E(int np, Database* db, int initial, int target, bool useFile);
void outputProtocol(int N, double T, int targetLump, const Eigen::VectorXd& init, 
										const Eigen::MatrixXd& probs, const Eigen::MatrixXd& prot, 
										std::ostream& ofile);






void testEqDeriv(int N, Database* db);
void testProbGradient(int np, Database* db);



}