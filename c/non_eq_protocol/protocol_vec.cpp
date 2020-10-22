

#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include "pair.h"
#include "adjacency.h"
#include "database.h"
#include "nauty.h"
#include "graph.h"
#include "graphviz.h"
#include "protocol.h"
#include "design.h"
#include "../defines.h"
namespace bd{

#define LOWER 0.1
#define UPPER 10000.0 
#define GAM exp(5.58)


/**********************************************************************/
/******************** Vector Valued Protocols *************************/
/**********************************************************************/



OptInfo::OptInfo(int N_, int num_states_, int numTypes_, int numInteractions_, int* particleTypes_, 
								 int initial_, int target_, double rho_, Database* db_) {
	problem_type = 0;
	N = N_; //here N is the number of particles, not temp discretization
	num_states = num_states_; initial = initial_; target = target_;
	numTypes = numTypes_; numInteractions = numInteractions_;
	rho = rho_; E = 0;
	h1 = 0; h2 = 0;
	db = db_;

	particleTypes = new int[N];
	for (int i = 0; i < N; i++) {
		particleTypes[i] = particleTypes_[i];
	}

}



void eigToDlib(int N, int numInteractions, const Eigen::MatrixXd& kappaMat, 
							column_vector& kappaVec) {
	//convert an eigen matrix of kappa values to a vector
	//matrix is num by N. vector has each interaction type stacked.

	int count = 0;
	for (int i = 0; i < numInteractions; i++) {
		for (int j = 0; j < N; j++) {
			kappaVec(count) = kappaMat(i,j);
			count++;
		}
	}

}

void dlibToEig(int N, int numInteractions, const column_vector kappaVec, 
							 Eigen::MatrixXd& kappaMat) {
	//convert from vector to matrix

	int count = 0;
	for (int i = 0; i < numInteractions; i++) {
		for (int j = 0; j < N; j++) {
			kappaMat(i,j) = kappaVec(count);
			count++;
		}
	}
}

void makeKappaMap(int numTypes, const Eigen::VectorXd& kappaVals, 
									std::map<std::pair<int,int>,double>& kappa) {
	//fill the map with values in kappaVals

	int counter = 0; 

	for (int i = 0; i < numTypes; i++) {
		for (int j = i; j < numTypes; j++) {
			kappa[std::make_pair(i,j)] = kappaVals(counter);
			kappa[std::make_pair(j,i)] = kappaVals(counter);
			counter++;
		}
	}
}

void applyBounds(int elements, Eigen::VectorXd& kappaVec, 
								 double lower, double upper) {
	//apply upper and lower bounds to each component of the vector

	for (int i = 0; i < elements; i++) {
		if (kappaVec(i) > upper) {
			kappaVec(i) = upper;
		}
		else if (kappaVec(i) < lower) {
			kappaVec(i) = lower;
		}
	}
}








void OptInfo::CreateRateMatrix(Eigen::MatrixXd& R, const Eigen::MatrixXd& r_form, 
															 const Eigen::VectorXd& kappaVals) {
	//use the database to construct the rate matrix at the given values
	//no derivative, so this is just a matrix

	//declare equilibrium measure and perform reweight
	double* eq = new double[num_states];           //equilibrium measure
	std::map<std::pair<int,int>,double> kappa;    //map for interaction type to kappa
	makeKappaMap(numTypes, kappaVals, kappa);
	reweight(N, num_states, db, particleTypes, eq, kappa);

	//fill in forward and satisfy detailed balance wrt eq for backward
	for (int i = 0; i < num_states; i++) {
		for (int j = 0; j < num_states; j++) {
			if (r_form(i,j) > 0) {
				R(i,j) = r_form(i,j);
				R(j,i) = r_form(i,j) * (eq[i] / eq[j]);
			}
		}
	}

	//finally, must set the diagonal to be negative of the row sum
	for (int i = 0; i < num_states; i++) {
		double R_sum = 0; 

		for (int j = 0; j < num_states; j++) {
			R_sum += R(i,j);
		}

		R(i,i) = -R_sum; 
	}


	delete []eq;
}

void eqDerivative(int N, int num_states, int numTypes, int numInteractions, Database* db, 
									int bond, const double* eq, double* eqDeriv, int* particleTypes, 
									const Eigen::VectorXd& kappaVals) {
	//compute the derivative of the eq measure wrt to the given kappa

	//map bond type to index
	std::map<std::pair<int,int>,int> index;
	int* indexVals = new int[numInteractions];
	for (int i = 0; i < numInteractions; i++) {
		indexVals[i] = i;
	}
	makeIndexMap(numTypes, indexVals, index);

	//init the weighted sum
	double S = 0; 

	//sum the eq prob over all states, weighted by n_b
	for (int i = 0; i < num_states; i++) {
		int b = searchForBonds(N, db, i, bond, particleTypes, index);
		S += eq[i] * float(b);
	}

	//fill in the derivative array
	for (int i = 0; i < num_states; i++) {
		int b = searchForBonds(N, db, i, bond, particleTypes, index);
		eqDeriv[i] = eq[i] / kappaVals[bond] * (float(b) - S);
	}

	delete []indexVals;

}

void OptInfo::CreateRateMatrix(Eigen::MatrixXd& R, Eigen::MatrixXd& D, 
															 const Eigen::MatrixXd& r_form, const Eigen::VectorXd& kappaVals, 
															 int ixn) {
	//use the database to construct the rate matrix and its kappa_i derivative

	//declare equilibrium measure and perform reweight
	double* eq = new double[num_states];           //equilibrium measure
	double* eqDeriv = new double[num_states];      //eq measure derivatives
	std::map<std::pair<int,int>,double> kappaMap;    //map for interaction type to kappa

	//make kappa mapping
	makeKappaMap(numTypes, kappaVals, kappaMap);

	//perform reqeighting of the equilirium measure
	reweight(N, num_states, db, particleTypes, eq, kappaMap);

	//compute the ixn gradient of eq measure
	eqDerivative(N, num_states, numTypes, numInteractions, db, ixn, eq, 
							 eqDeriv, particleTypes, kappaVals);


	//fill in forward and satisfy detailed balance wrt eq for backward
	for (int i = 0; i < num_states; i++) {
		for (int j = 0; j < num_states; j++) {
			if (r_form(i,j) > 0) {
				//forward entries
				R(i,j) = r_form(i,j); D(i,j) = 0;
				//backward entries
				R(j,i) = r_form(i,j) * (eq[i] / eq[j]);
				D(j,i) = r_form(i,j) * (eqDeriv[i] * eq[j] - eq[i] * eqDeriv[j]) / (eq[j]*eq[j]);
			}
		}
	}

	//finally, must set the diagonal to be negative of the row sum
	for (int i = 0; i < num_states; i++) {
		double R_sum = 0; double D_sum = 0;

		for (int j = 0; j < num_states; j++) {
			R_sum += R(i,j); D_sum += D(i,j);
		}

		R(i,i) = -R_sum; D(i,i) = -D_sum;
	}


	delete []eq; delete []eqDeriv;
	
}

void evalGradient(int N, const Eigen::MatrixXd& kappaMat, int num_states, 
							 OptInfo* problem, const Eigen::MatrixXd& f_rate, const Eigen::VectorXd& t_disc, 
							 const Eigen::VectorXd& t_centers, int initial, std::vector<int> targets,
							 Eigen::VectorXd& R) {

	//loop over each interaction type
	int numInteractions = problem->GetNumInteractions();
	for (int ixn = 0; ixn < numInteractions; ixn++) {

		//loop over each component of R, compute the gradient
		for (int j = 0; j < N; j++) {
			//compute the derivative of the rate matrix at current beta val
			Eigen::MatrixXd rate; rate.setZero(num_states,num_states);
			Eigen::MatrixXd LD; LD.setZero(num_states,num_states);
			problem->CreateRateMatrix(rate, LD, f_rate, kappaMat.col(j), ixn);
			//std::cout << rate << "\n";
			//std::cout << rate(7,7) << "\n";
			//Eigen::MatrixXd test(num_states,num_states); test = rate * 100; test = test.exp();
			//std::cout << test << "\n";
			//abort();

			//apply the transition operator from 0 to t_center(j), then LD, then to T
			Eigen::MatrixXd T1; T1.setZero(num_states,num_states);
			Eigen::MatrixXd T2; T2.setZero(num_states,num_states);
			Eigen::MatrixXd trans_op; trans_op.setZero(num_states,num_states);
			
			createTransitionOperator(N, num_states, problem, f_rate, kappaMat, t_disc, 0, 
															t_centers(j), j, T1);
			createTransitionOperator(N, num_states, problem, f_rate, kappaMat, t_disc, t_centers(j), 
															t_disc(N), j, T2);
			trans_op = T1 * LD * T2;
			//std::cout << LD.row(initial) << "\n";
			//abort();

			//std::cout << trans_op << "\n";


			//extract the relevant entry: row = initial, column = target
			for (int t = 0; t < targets.size(); t++) {
				R(ixn*N+j) += trans_op(initial, targets[t]);
			}

			//check for a move that violates the lower bound
			//todo - not neccesary with dlib constraint opt?
		}
	}

}

void evalGradientFD(int N, Eigen::MatrixXd& kappaMat, int num_states, 
							 OptInfo* problem, const Eigen::MatrixXd& f_rate, const Eigen::VectorXd& t_disc, 
							 const Eigen::VectorXd& t_centers, int initial, std::vector<int> targets,
							 Eigen::VectorXd& R) {

	//evaluate the gradient by finite differencing

	double h = 1e-6;
	//get the current probability once
	double p0 = computeFinalProb(N, num_states, problem, f_rate, kappaMat, t_disc,
														  initial, targets);

	//loop over every element of kappaMat and compute the partial deriv by fd

	//loop over each interaction type
	int numInteractions = problem->GetNumInteractions();
	for (int ixn = 0; ixn < numInteractions; ixn++) {

		//loop over each component of R, compute the gradient
		for (int j = 0; j < N; j++) {
			
			//increment relevant component by h
			kappaMat(ixn, j) += h;
			//get the forward value
			double p1 = computeFinalProb(N, num_states, problem, f_rate, kappaMat, t_disc,
														  initial, targets);
			//go back to original mat
			kappaMat(ixn, j) -= h;

			//compute finite difference
			R(ixn*N+j) = (p1-p0) / h;

		}
	}

}

Eigen::RowVectorXd push_forward(int num_states, double dt, Eigen::MatrixXd& rate,  
									const Eigen::RowVectorXd& initial) {
	//apply the given rate matrix forward for the given time step

	rate *= dt;
	//exponentiate to get transition operator
	Eigen::MatrixXd trans_op; 
	trans_op.setZero(num_states,num_states);
	trans_op = rate.exp();

	//do the matrix vector mult
	Eigen::RowVectorXd final; final.setZero(num_states);
	final = initial * trans_op;

	return final;
}

Eigen::VectorXd push_back(int num_states, double dt, Eigen::MatrixXd& rate,  
									const Eigen::VectorXd& initial) {
	//apply the given rate matrix forward for the given time step

	rate *= dt;
	//exponentiate to get transition operator
	Eigen::MatrixXd trans_op; 
	trans_op.setZero(num_states,num_states);
	trans_op = rate.exp();

	//do the matrix vector mult
	Eigen::VectorXd final; final.setZero(num_states);
	final = trans_op * initial;

	return final;
}

void f_solve(int N, Eigen::MatrixXd& kappaMat, int num_states, OptInfo* problem, 
						 const Eigen::MatrixXd& f_rate, const Eigen::VectorXd& t_disc, 
						 Eigen::MatrixXd& f_solution) {
	//solve forward kolm. equation

	//init the rate matrix 
	Eigen::MatrixXd rate; 

	for (int i = 0; i < N; i++) {
		//compute the rate matrix at this interval
		rate.setZero(num_states,num_states);
		problem->CreateRateMatrix(rate, f_rate, kappaMat.col(i));
		//multiply matrix by timestep
		double delta_t = t_disc(i+1) - t_disc(i);

		//push_forward to the next time
		f_solution.row(i+1) = push_forward(num_states, delta_t, rate, f_solution.row(i));
	}


}

void f_solve_half(int N, Eigen::MatrixXd& kappaMat, int num_states, OptInfo* problem, 
						 const Eigen::MatrixXd& f_rate, const Eigen::VectorXd& t_disc, 
						 Eigen::MatrixXd& f_solution) {
	//solve forward kolm. equation

	//init the rate matrix 
	Eigen::MatrixXd rate; 

	for (int i = 0; i < N; i++) {
		//compute the rate matrix at this interval
		rate.setZero(num_states,num_states);
		problem->CreateRateMatrix(rate, f_rate, kappaMat.col(i));
		//multiply matrix by timestep
		double delta_t = t_disc(i+1) - t_disc(i);
		if (i == 0) delta_t /= 2.0;

		//multiply previous probability vector by transition operator
		f_solution.row(i+1) = push_forward(num_states, delta_t, rate, f_solution.row(i));
	}


}

void b_solve(int N, Eigen::MatrixXd& kappaMat, int num_states, OptInfo* problem, 
						 const Eigen::MatrixXd& f_rate, const Eigen::VectorXd& t_disc, 
						 Eigen::MatrixXd& b_solution) {
	//solve forward kolm. equation

	//init the rate matrix 
	Eigen::MatrixXd rate; 

	for (int i = N; i > 0; i--) {
		//compute the rate matrix at this interval
		rate.setZero(num_states,num_states);
		problem->CreateRateMatrix(rate, f_rate, kappaMat.col(i-1));
		//multiply matrix by timestep
		double delta_t = t_disc(i) - t_disc(i-1);

		//multiply previous probability vector by transition operator
		b_solution.col(i-1) = push_back(num_states, delta_t, rate, b_solution.col(i));
	}


}

void b_solve_half(int N, Eigen::MatrixXd& kappaMat, int num_states, OptInfo* problem, 
						 const Eigen::MatrixXd& f_rate, const Eigen::VectorXd& t_disc, 
						 Eigen::MatrixXd& b_solution) {
	//solve forward kolm. equation

	//init the rate matrix 
	Eigen::MatrixXd rate; 

	for (int i = N; i > 0; i--) {
		//compute the rate matrix at this interval
		rate.setZero(num_states,num_states);
		problem->CreateRateMatrix(rate, f_rate, kappaMat.col(i-1));
		//multiply matrix by timestep
		double delta_t = t_disc(i) - t_disc(i-1);
		if (i == N) delta_t /= 2.0;
		
		//multiply previous probability vector by transition operator
		b_solution.col(i-1) = push_back(num_states, delta_t, rate, b_solution.col(i));
	}


}

double gaussQuad(double t0, double tf, int num_states, const Eigen::MatrixXd& rate,
								 const Eigen::MatrixXd& LD, const Eigen::RowVectorXd& pL,
								 const Eigen::VectorXd muR) {
	//use gaussian quadrature to evaluate the integrals

	Eigen::MatrixXd baseRate; baseRate.setZero(num_states,num_states);

	//midpoint and interval length
	double mid = (t0+tf)/2.0;
	double dt = tf-t0;

	/*
	//define the points
	double x0 = mid - (0.5*dt) / sqrt(3); double x1 = mid + (0.5*dt) / sqrt(3);

	//get p at quadrature points
	Eigen::RowVectorXd p0; p0.setZero(num_states); 
	Eigen::RowVectorXd p1; p1.setZero(num_states);

	baseRate = rate; p0 = push_forward(num_states, x0-t0, baseRate, pL);
	baseRate = rate; p1 = push_forward(num_states, x1-t0, baseRate, pL);

	//get mu at quadrature points
	Eigen::VectorXd mu0; mu0.setZero(num_states); 
	Eigen::VectorXd mu1; mu1.setZero(num_states);

	baseRate = rate; mu0 = push_back(num_states, tf-x0, baseRate, muR);
	baseRate = rate; mu1 = push_back(num_states, tf-x1, baseRate, muR);

	double I_0 = p0 * LD * mu0;
	double I_1 = p1 * LD * mu1;
	return 0.5*dt*(I_0+I_1);
	*/

	/*
	//define the points and weights
	double xm = mid - (0.5*dt)*sqrt(3.0/5.0); double wm = 5.0/9.0;
	double xp = mid + (0.5*dt)*sqrt(3.0/5.0); double wp = 5.0/9.0;
	double x0 = mid; double w0 = 8.0/9.0;

	//get p at quadrature points
	Eigen::RowVectorXd pm; pm.setZero(num_states); 
	Eigen::RowVectorXd pp; pp.setZero(num_states);
	Eigen::RowVectorXd p0; p0.setZero(num_states);

	baseRate = rate; pm = push_forward(num_states, xm-t0, baseRate, pL);
	baseRate = rate; pp = push_forward(num_states, xp-t0, baseRate, pL);
	baseRate = rate; p0 = push_forward(num_states, x0-t0, baseRate, pL);

	//get mu at quadrature points
	Eigen::VectorXd mum; mum.setZero(num_states); 
	Eigen::VectorXd mup; mup.setZero(num_states);
	Eigen::VectorXd mu0; mu0.setZero(num_states);

	baseRate = rate; mum = push_back(num_states, tf-xm, baseRate, muR);
	baseRate = rate; mup = push_back(num_states, tf-xp, baseRate, muR);
	baseRate = rate; mu0 = push_back(num_states, tf-x0, baseRate, muR);

	//std::cout << muR << "\n";
	//std::cout << mum << "\n";

	double I_m = pm * LD * mum;
	double I_p = pp * LD * mup;
	double I_0 = p0 * LD * mu0;
	return 0.5*dt*(I_m * wm + I_p * wp + I_0 * w0);
	*/

	int G = 8;
	Eigen::VectorXd x; x.setZero(G); Eigen::VectorXd w; w.setZero(G);
	x(0) = -0.9602; x(1) = -0.7966; x(2) = -0.5255; x(3) = -0.1834; 
	x(7) = -x(0); x(6) = -x(1); x(5) = -x(2); x(4) = -x(3);

	w(0) = 0.101228; w(1) = 0.22238; w(2) = 0.3137; w(3) = 0.36268;
	w(7) = w(0); w(6) = w(1); w(5) = w(2); w(4) = w(3);

	double I = 0;

	for (int i = 0; i < G; i++) {
		double final_time = mid + 0.5*dt*x(i);

		Eigen::RowVectorXd p; p.setZero(num_states); 
		baseRate = rate; p = push_forward(num_states, final_time-t0, baseRate, pL);

		Eigen::VectorXd mu; mu.setZero(num_states); 
		baseRate = rate; mu = push_back(num_states, tf-final_time, baseRate, muR);

		double inc = p * LD * mu;
		I += w(i) * inc;
	}

	return 0.5 * dt * I;




}



void evalGradientAdjoint(int N, Eigen::MatrixXd& kappaMat, int num_states, 
							 OptInfo* problem, const Eigen::MatrixXd& f_rate, const Eigen::VectorXd& t_disc, 
							 const Eigen::VectorXd& t_centers, int initial, std::vector<int> targets,
							 Eigen::VectorXd& R) {
	//evaluate the gradient by using the adjoint formulation

	//first, we need to solve the forward and backward kolmogorov equations

	//init storage for solution vectors at each time
	Eigen::MatrixXd f_solution; f_solution.setZero(N+1, num_states); //rows are probabilities
	Eigen::MatrixXd f_solution_half; f_solution_half.setZero(N+1, num_states); //rows are probabilities
	Eigen::MatrixXd b_solution; b_solution.setZero(num_states, N+1); //columns are expectations
	Eigen::MatrixXd b_solution_half; b_solution_half.setZero(num_states, N+1); //columns are expectations

	//set the initial condition for forward solve
	f_solution(0, initial) = 1.0; f_solution_half(0, initial) = 1.0;

	//set the terminal condition for backward solve
	for (int i = 0; i < targets.size(); i++) {
		b_solution(targets[i], N) = 1.0;
		b_solution_half(targets[i], N) = 1.0;
	}

	//do the forward solve
	f_solve(N, kappaMat, num_states, problem, f_rate, t_disc, f_solution);
	f_solve_half(N, kappaMat, num_states, problem, f_rate, t_disc, f_solution_half);

	//do the backward solve
	b_solve(N, kappaMat, num_states, problem, f_rate, t_disc, b_solution);
	b_solve_half(N, kappaMat, num_states, problem, f_rate, t_disc, b_solution_half);

	
	//loop over each interaction, computing derivatives of the rate matrix, to get gradient
	int numInteractions = problem->GetNumInteractions();
	for (int ixn = 0; ixn < numInteractions; ixn++) {

		//loop over each component of R, compute the gradient
		for (int j = 0; j < N; j++) {
			//create the rate matrix and its derivative at this time
			Eigen::MatrixXd rate; rate.setZero(num_states,num_states);
			Eigen::MatrixXd LD; LD.setZero(num_states,num_states);
			problem->CreateRateMatrix(rate, LD, f_rate, kappaMat.col(j), ixn);

			//compute adjoint form of the gradient
			double g = (f_solution.row(j) * LD * b_solution.col(j)); //using left endpoints
			double gh = (f_solution_half.row(j+1) * LD * b_solution_half.col(j));
			double g2 = (f_solution.row(j+1) * LD * b_solution.col(j+1));
			//double g = f_solution.row(j+1) * LD * b_solution.col(j); //using midpoints

			//store the gradient at this time
			double dt = t_disc(j+1) - t_disc(j);
			double val = (g+4*gh+g2)/6.0 * dt;
			if (j == 0 || j == N-1) { //integral approx is bad at endpoints, use FD
				//val = gaussQuad(t_disc(j), t_disc(j+1), num_states, rate ,LD,
														  //f_solution.row(j), b_solution.col(j+1));
				double h = 1e-7;
				double base = computeFinalProb(N, num_states, problem, f_rate, kappaMat, t_disc,
														  initial, targets);
				kappaMat(ixn,j) += h;
				double pert = computeFinalProb(N, num_states, problem, f_rate, kappaMat, t_disc,
														  initial, targets);
				kappaMat(ixn,j) -= h;
				val = (pert-base) / h;
			}
			//double val = (g+g2)/2.0 * dt;
			R(ixn*N+j) = val;

			//if the gradient would push protocols past bounds, set to 0
			if ((kappaMat(ixn,j) == UPPER && val > 0) ||
				  (kappaMat(ixn,j) == LOWER && val < 0)) {
				R(ixn*N+j) = 0;
			}

			//printf("g %f, gh %f, g2 %f, grad %f\n", g, gh, g2, R(ixn*N+j));
		}

	}

}

void evalGradientAdjoint(int N, Eigen::VectorXd& kappaVec, int num_states, 
							 OptInfo* problem, const Eigen::MatrixXd& f_rate, const Eigen::VectorXd& t_disc, 
							 const Eigen::VectorXd& t_centers, int initial, std::vector<int> targets,
							 Eigen::VectorXd& R) {
	//convert the vector to matrix to use the other overloaded function

	int numInteractions = problem->GetNumInteractions();
	Eigen::MatrixXd kappaMat; kappaMat.setZero(numInteractions, N);

	for (int ixn = 0; ixn < numInteractions; ixn++) {
		for (int j = 0; j < N; j++) {
			kappaMat(ixn, j) = kappaVec(ixn*N+j);
		}
	}

	evalGradientAdjoint(N, kappaMat, num_states, problem, f_rate, t_disc, t_centers, 
											initial, targets, R);
}


double protocolOptLineSearch(int N, Eigen::VectorXd& kappaVec, int num_states, 
							 OptInfo* problem, const Eigen::MatrixXd& f_rate, const Eigen::VectorXd& t_disc, 
							 const Eigen::VectorXd& t_centers, int initial, std::vector<int> targets,
							 Eigen::VectorXd& R, double p0) {
	//perform line search to determine a step size that guarantees progress
	//return the step size

	//parameters to the line search
	double h = 100.0;
	double alpha = 0.5;
	double attempts = 14;
	double tol = 1e-12;

	//parameters for new step test vector
	int numInteractions = problem->GetNumInteractions();
	Eigen::VectorXd testR;

	for (int i = 0; i < attempts; i++) {
		//at each attempt, check if step size improves objective. if not, reduce by alpha

		testR.setZero(N*numInteractions);
		testR = kappaVec + h * R;
		applyBounds(N*numInteractions, testR, LOWER, UPPER);
		double p1 = computeFinalProb(N, num_states, problem, f_rate, testR, t_disc,  
																initial, targets);
		if ((p1 - p0) > tol) {
			printf("The step has been reduced %d times\n", i);
			return h;
		} 
		else {
			h *= alpha;
		}
	}

	//if we get to this point, we did not see an increase in objective. return smallest h
	return h;

}

double protocolOptLineSearch_energy(int N, Eigen::VectorXd& energyVec, int num_states, 
							 OptInfo* problem, const Eigen::MatrixXd& f_rate, const Eigen::VectorXd& t_disc, 
							 const Eigen::VectorXd& t_centers, int initial, std::vector<int> targets,
							 Eigen::VectorXd& R, double p0) {
	//perform line search to determine a step size that guarantees progress
	//return the step size

	//parameters to the line search
	double h = 5.0;
	double alpha = 0.5;
	double attempts = 14;
	double tol = 1e-12;

	//parameters for new step test vector
	int numInteractions = problem->GetNumInteractions();
	Eigen::VectorXd testR;

	for (int i = 0; i < attempts; i++) {
		//at each attempt, check if step size improves objective. if not, reduce by alpha

		testR.setZero(N*numInteractions);
		testR = energyVec + h * R;
		applyBounds(N*numInteractions, testR, 1.0, 14.0);
		for (int j = 0; j < N*numInteractions; j++) {
			testR(j) = exp(testR(j)) / sqrt(2.0*40*40*testR(j));
		}
		double p1 = computeFinalProb(N, num_states, problem, f_rate, testR, t_disc,  
																initial, targets);
		//printf("Attempt %i, prob %f\n", i, p1);
		if ((p1 - p0) > tol) {
			printf("The step has been reduced %d times\n", i);
			return h;
		} 
		else {
			h *= alpha;
		}
	}

	//if we get to this point, we did not see an increase in objective. return smallest h
	return h;

}

void diffusionSolve(int N, int numInteractions, double lambda, double h, double dt,
										const Eigen::VectorXd& pre, Eigen::VectorXd& protocol) {
	//solve the diffusion equation corresponding to regularizing term
	//neumann boundary conditions


	//define the diffusive cfl
	h = dt;
  double cfl = lambda * h / (dt*dt);      

	//loop over each interaction, solving the diffusion equation
	for (int ixn = 0; ixn < numInteractions; ixn++) {

		//rhs and solution vectors
		Eigen::VectorXd rhs; rhs.setZero(N); rhs = protocol.segment(ixn*N,N);
		Eigen::VectorXd mu; mu.setZero(N);   mu  = pre.segment(ixn*N,N) * cfl;
		Eigen::VectorXd soln; soln.setZero(N);

		//storage for sweep vectors for Thomas algo
		Eigen::VectorXd c; c.setZero(N-1); c(0) = - 2.0 * mu(0) / (1.0+2.0*mu(0));
		Eigen::VectorXd d; d.setZero(N);   d(0) = rhs(0) / (1.0+2.0*mu(0));

		//apply the Thomas algo to solve the tridiagonal system
		for (int i = 1; i < N-1; i++) {
			double denom = 1.0 + 2.0*mu(i) + mu(i)*c(i-1);
			c(i) = -mu(i) / denom;
			d(i) = (rhs(i) + mu(i)*d(i-1)) / denom;
		}
		d(N-1) = (rhs(N-1) + 2.0*mu(N-1)*d(N-2)) / (1.0+2.0*mu(N-1) +2.0*mu(N-1)*c(N-2));

		soln(N-1) = d(N-1);
		for (int i = N-2; i > -1; i--) {
			soln(i) = d(i) - c(i) * soln(i+1);
		} 

		//put soln back in protocol vector
		protocol.segment(ixn*N,N) = soln;
	}
}


void protocolOptAdjoint(int N, Eigen::MatrixXd& kappaMat, int num_states, 
							 OptInfo* problem, const Eigen::MatrixXd& f_rate, const Eigen::VectorXd& t_disc, 
							 const Eigen::VectorXd& t_centers, int initial, std::vector<int> targets) {
	//iterative algorithm to determine optimal protocol using adjoint method

	//gradient descent parameters
	int max_iter = 40;
	bool regular = false;
	double dt = t_disc(1) - t_disc(0);
	double T = t_disc(N);
	double lambda = dt;

	//init vector to store gradient
	int numInteractions = problem->GetNumInteractions();
	Eigen::VectorXd R; 

	//convert kappaMat to vector storage
	Eigen::VectorXd kappaVec; kappaVec.setZero(N*numInteractions);
	for (int ixn = 0; ixn < numInteractions; ixn++) {
		kappaVec.segment(ixn*N, N) = kappaMat.row(ixn).transpose();
	} 

	//get the objective value at the IC
	double p0 = computeFinalProb(N, num_states, problem, f_rate, kappaVec, t_disc,  
																initial, targets);

	//do the iteration
	for (int iter = 0; iter < max_iter; iter++ ) {
		//first, compute the gradient using the adjoint
		R.setZero(N*numInteractions);
		evalGradientAdjoint(N, kappaVec, num_states, problem, f_rate, t_disc, t_centers, 
											initial, targets, R);

		//normalize and do line search for step size
		double r_norm = R.norm();
		R /= r_norm;
		double h = protocolOptLineSearch(N, kappaVec, num_states, problem, f_rate, t_disc, 
							 											 t_centers, initial, targets, R, p0);

		//take the gradient descent step
		kappaVec += h * R;
		applyBounds(N*numInteractions, kappaVec, LOWER, UPPER);
		std::cout << R << "\n";

		//printf("Before the diffusion solve we have\n");
		//std::cout << kappaVec << "\n";

		//check if pre-conditioner should be defined
		bool precondition = false;
		Eigen::VectorXd pre; pre.setOnes(N*numInteractions);
		if (precondition) {
			//make function to compute diagonal preconditioner

		}

		//add in the regularizer if desired
		if (regular) {
			//add in the regulariztion term
			diffusionSolve(N, numInteractions, h, lambda, dt, pre, kappaVec);
		}

		//printf("After the diffusion solve we have\n");
		//std::cout << kappaVec << "\n";

		//evaluate the objective at the new values, set to p0
		p0 = computeFinalProb(N, num_states, problem, f_rate, kappaVec, t_disc,  
																initial, targets);

		//print out the progress
		printf("Iteration: %d, Objective Function Value: %f, Gradient Norm: %f\n", iter, p0, r_norm);


		


	}

	//return from vector format back to matrix
	for (int ixn = 0; ixn < numInteractions; ixn++) {
		kappaMat.row(ixn) = kappaVec.segment(ixn*N,N).transpose();
	}

	std::cout << kappaMat << "\n";




}

void protocolOptAdjoint_energy(int N, Eigen::MatrixXd& kappaMat, int num_states, 
							 OptInfo* problem, const Eigen::MatrixXd& f_rate, const Eigen::VectorXd& t_disc, 
							 const Eigen::VectorXd& t_centers, int initial, std::vector<int> targets) {
	//iterative algorithm to determine optimal protocol using adjoint method

	//algorithm parameters
	bool regular = true;
	bool line_search = false;
	bool precondition = true;

	//gradient descent parameters
	int max_iter = 100;
	double dt = t_disc(1) - t_disc(0);
	double T = t_disc(N);
	double lambda = dt;

	//init vector to store gradient
	int numInteractions = problem->GetNumInteractions();
	Eigen::VectorXd R; 

	//convert kappaMat to vector storage
	Eigen::VectorXd kappaVec; kappaVec.setZero(N*numInteractions);
	for (int ixn = 0; ixn < numInteractions; ixn++) {
		kappaVec.segment(ixn*N, N) = kappaMat.row(ixn).transpose();
	} 

	//make a vector of energies as well
	Eigen::VectorXd energyVec; energyVec.setZero(N*numInteractions);
	for (int i = 0; i < N*numInteractions; i++) {
		energyVec(i) = stickyNewton(8.0, 40.0, kappaVec(i), 1.0); 
	}

	//get the objective value at the IC
	double p0 = computeFinalProb(N, num_states, problem, f_rate, kappaVec, t_disc,  
																initial, targets);

	//do the iteration
	for (int iter = 0; iter < max_iter; iter++ ) {
		//first, compute the gradient using the adjoint
		R.setZero(N*numInteractions);
		evalGradientAdjoint(N, kappaVec, num_states, problem, f_rate, t_disc, t_centers, 
											initial, targets, R);

		//multiply by d(kappa)/d(E)
		for (int element = 0; element < N*numInteractions; element++) {
			double E = energyVec(element);
			double chain_rule = exp(E) / (2.0*40.0*E) * ( (2.0*E-1.0) / sqrt(2.0*E));
			R(element) *= chain_rule;
		}

		//check if line search should be done
		double r_norm = R.norm();
		double h;
		if (line_search) {
			R /= r_norm;
			h = protocolOptLineSearch_energy(N, energyVec, num_states, problem, f_rate, t_disc, 
							 											 t_centers, initial, targets, R, p0);
		}
		else {
			h = dt;
			R /= h;
		}

		//check if pre-conditioner should be defined
		Eigen::VectorXd pre; pre.setOnes(N*numInteractions);
		if (precondition) {
			//make function to compute diagonal preconditioner

			//try something
			for (int i = 0; i < N*numInteractions; i++) {
				pre(i) = (N-(i%N))*10;
			}
		}

		//take the gradient descent step
		energyVec += h * pre.cwiseProduct(R);
		applyBounds(N*numInteractions, energyVec, 1.0, 14.0);
		//std::cout << R << "\n";

		//add in the regularizer if desired
		if (regular) {
			//add in the regulariztion term
			diffusionSolve(N, numInteractions, h, lambda, dt, pre, energyVec);
		}

		//evaluate the objective at the new values, set to p0
		for (int i = 0; i < N*numInteractions; i++) {
			kappaVec(i) = exp(energyVec(i)) / sqrt(2.0*40*40*energyVec(i));
		}
		p0 = computeFinalProb(N, num_states, problem, f_rate, kappaVec, t_disc,  
																initial, targets);

		//print out the progress
		printf("Iteration: %d, Objective Function Value: %f, Gradient Norm: %f\n", iter, p0, r_norm);


		


	}

	//return from vector format back to matrix
	for (int ixn = 0; ixn < numInteractions; ixn++) {
		kappaMat.row(ixn) = kappaVec.segment(ixn*N,N).transpose();
	}

	std::cout << kappaMat << "\n";




}









void createTransitionOperator(int N, int num_states, OptInfo* problem, 
															const Eigen::MatrixXd& f_rate, const Eigen::MatrixXd& kappaMat,
															const Eigen::VectorXd& t_disc, double t0, double tf, int interval, 
															Eigen::MatrixXd& trans_op) {
	//construct the transition operator from t0 to tf, given the protocol

	//begin by setting the transition operator equal to the identity
	trans_op = Eigen::MatrixXd::Identity(num_states,num_states);

	//define a storage for the rate matrices
	Eigen::MatrixXd rate; 

	//check if t0 is 0 and tf is not T. if yes, propogate until tf
	if (t0 == 0 && tf != t_disc(N)) {
		//the final index is interval
		int final = interval; 
		//printf("Final interval is %d\n", interval);
		//printf("Flag1\n");

		//loop over each interval that gets the full time step
		for (int i = 0; i < final; i++) {
			rate.setZero(num_states,num_states);
			problem->CreateRateMatrix(rate, f_rate, kappaMat.col(i));;
			double delta_t = t_disc(i+1) - t_disc(i);
			rate *= delta_t;
			trans_op *= rate.exp();
		}

		//the final interval gets a partial step
		rate.setZero(num_states,num_states);
		problem->CreateRateMatrix(rate, f_rate, kappaMat.col(final));
		double delta_t = tf - t_disc(final);
		rate *= delta_t;
		trans_op *= rate.exp();

		return;
	}

	//check if tf is the final time and t0 is not 0
	if (tf == t_disc(N) && t0 != 0) {
		//int initial index is interval
		int first = interval;
		//printf("Flag2\n");
		//printf("First interval is %d\n", interval);

		//the first interval gets a partial step
		rate.setZero(num_states,num_states);
		problem->CreateRateMatrix(rate, f_rate, kappaMat.col(first));
		double delta_t = t0 - t_disc(first);
		rate *= delta_t;
		trans_op *= rate.exp();

		//loop over the rest of the intervals
		for (int i = first+1; i < N; i++) {
			rate.setZero(num_states,num_states);
			problem->CreateRateMatrix(rate, f_rate, kappaMat.col(i));
			double delta_t = t_disc(i+1) - t_disc(i);
			rate *= delta_t;
			trans_op *= rate.exp();
		}

		return;
	}

	//check if t0 == 0 and tf == T
	if (tf == t_disc(N) && t0 == 0) {

		//printf("Flag3\n");
		for (int i = 0; i < N; i++) {
			rate.setZero(num_states,num_states);
			problem->CreateRateMatrix(rate, f_rate, kappaMat.col(i));
			double delta_t = t_disc(i+1) - t_disc(i);
			rate *= delta_t;
			trans_op *= rate.exp();
		}

		return;
	}

}





double computeFinalProb(int N, int num_states, OptInfo* problem, const Eigen::MatrixXd& f_rate,
												const Eigen::MatrixXd& kappaMat, const Eigen::VectorXd& t_disc,  
												int initial, std::vector<int> targets) {

	Eigen::MatrixXd trans_op = Eigen::MatrixXd::Identity(num_states,num_states);
	Eigen::MatrixXd rate; rate.setZero(num_states,num_states);

	for (int i = 0; i < N; i++) {
		//create the rate matrix
		rate.setZero(num_states,num_states);
		problem->CreateRateMatrix(rate, f_rate, kappaMat.col(i));

		//scale by time interval and exponentiate
		double delta_t = t_disc(i+1) - t_disc(i);
		rate *= delta_t;
		trans_op *= rate.exp();
	}

	double p = 0;
	for (int i = 0; i < targets.size(); i++) {
		p += trans_op(initial, targets[i]);
	}

	return p;
}

double computeFinalProb(int N, int num_states, OptInfo* problem, const Eigen::MatrixXd& f_rate,
												const Eigen::VectorXd& kappaVec, const Eigen::VectorXd& t_disc,  
												int initial, std::vector<int> targets) {

	Eigen::MatrixXd trans_op = Eigen::MatrixXd::Identity(num_states,num_states);
	Eigen::MatrixXd rate; rate.setZero(num_states,num_states);
	int numInteractions = problem->GetNumInteractions();
	Eigen::VectorXd theta; theta.setZero(numInteractions);

	for (int i = 0; i < N; i++) {
		//create the rate matrix
		rate.setZero(num_states,num_states);
		for (int ixn = 0; ixn < numInteractions; ixn++) {
			theta(ixn) = kappaVec(N*ixn+i);
		}
		problem->CreateRateMatrix(rate, f_rate, theta);

		//scale by time interval and exponentiate
		double delta_t = t_disc(i+1) - t_disc(i);
		rate *= delta_t;
		trans_op *= rate.exp();
	}

	double p = 0;
	for (int i = 0; i < targets.size(); i++) {
		p += trans_op(initial, targets[i]);
	}

	return p;
}

void computeProbVec(int N, int num_states, OptInfo* problem, const Eigen::MatrixXd& f_rate,
												const Eigen::MatrixXd& kappaMat, const Eigen::VectorXd& t_disc,  
												int initial, std::vector<int> targets, Eigen::MatrixXd& probs) {

	Eigen::MatrixXd trans_op = Eigen::MatrixXd::Identity(num_states,num_states);
	Eigen::MatrixXd rate; rate.setZero(num_states,num_states);

	for (int i = 0; i < N; i++) {
		rate.setZero(num_states,num_states);
		problem->CreateRateMatrix(rate, f_rate, kappaMat.col(i));

		//compute the eigenvalues of the rate matrix
		Eigen::EigenSolver<Eigen::MatrixXd> es(rate);
		//get second largest eigenvalue
		double e = es.eigenvalues()[0].real(); if (fabs(e) < 1e-8) e = es.eigenvalues()[1].real();
		for (int j = 1; j < num_states; j++) {
			double val = es.eigenvalues()[j].real();
			if (val > e && fabs(val) > 1e-8) {
				e = val;
			}
		}
		
		std::cout << "Eigenvalue 1 is:" << std::endl << e << std::endl;
		

		double delta_t = t_disc(i+1) - t_disc(i);
		rate *= delta_t;
		trans_op *= rate.exp();

		probs.row(i) = trans_op.row(initial);
	}

}



double dlib_objective_call(int N, int num_states, int numInteractions, OptInfo* problem, 
													 column_vector& kappaVec, const Eigen::MatrixXd& f_rate, 
													 const Eigen::VectorXd& t_disc, const Eigen::VectorXd& t_centers, 
													 int initial, std::vector<int> targets) {
	//call the dlib optimizer using an objective delta stop condition

	double m = dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),  // Use BFGS search algorithm
             dlib::objective_delta_stop_strategy(1e-14,40).be_verbose(), 

             [&](const column_vector& b) {
             	Eigen::MatrixXd kappaMat; kappaMat.setZero(numInteractions, N); 
             	dlibToEig(N, numInteractions, b, kappaMat);
             	return -computeFinalProb(N, num_states, problem, f_rate,
             	kappaMat, t_disc, initial, targets);}, 

             [&](const column_vector& b) {
             	//compute gradient 
             	Eigen::MatrixXd kappaMat; kappaMat.setZero(numInteractions, N); 
             	dlibToEig(N, numInteractions, b, kappaMat);
             	Eigen::VectorXd	R; R.setZero(N*numInteractions);
             	evalGradient(N, kappaMat, num_states, problem, f_rate, t_disc, t_centers, 
             							 initial, targets, R);
             	column_vector G(N*numInteractions); 
             	for (int i = 0; i < N*numInteractions; i++) G(i) = -R(i);
             	return G;
             }, 
             kappaVec, 0.1, 10000.0);

	return m;

}

double dlib_objective_call_fd(int N, int num_states, int numInteractions, OptInfo* problem, 
													 column_vector& kappaVec, const Eigen::MatrixXd& f_rate, 
													 const Eigen::VectorXd& t_disc, const Eigen::VectorXd& t_centers, 
													 int initial, std::vector<int> targets) {
	//call the dlib optimizer using an objective delta stop condition
	//finite difference for gradient

		double m = dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),  // Use BFGS search algorithm
             dlib::objective_delta_stop_strategy(1e-14,40).be_verbose(), 

             [&](const column_vector& b) {
             	Eigen::MatrixXd kappaMat; kappaMat.setZero(numInteractions, N); 
             	dlibToEig(N, numInteractions, b, kappaMat);
             	return -computeFinalProb(N, num_states, problem, f_rate,
             	kappaMat, t_disc, initial, targets);}, 

             [&](const column_vector& b) {
             	//compute gradient 
             	Eigen::MatrixXd kappaMat; kappaMat.setZero(numInteractions, N); 
             	dlibToEig(N, numInteractions, b, kappaMat);
             	Eigen::VectorXd	R; R.setZero(N*numInteractions);
             	evalGradientFD(N, kappaMat, num_states, problem, f_rate, t_disc, t_centers, 
             							 initial, targets, R);
             	column_vector G(N*numInteractions); 
             	for (int i = 0; i < N*numInteractions; i++) G(i) = -R(i);
             	return G;
             }, 
             kappaVec, 0.1, 10000.0);

	return m;

}

double dlib_objective_call_ng(int N, int num_states, int numInteractions, OptInfo* problem, 
													 column_vector& kappaVec, const Eigen::MatrixXd& f_rate, 
													 const Eigen::VectorXd& t_disc, const Eigen::VectorXd& t_centers, 
													 int initial, std::vector<int> targets) {
	//call the dlib optimizer using an objective delta stop condition
	//gradient free

	double m = dlib::find_min_using_approximate_derivatives(dlib::bfgs_search_strategy(),  // Use BFGS search algorithm
             dlib::objective_delta_stop_strategy(1e-14,40).be_verbose(), 

             [&](const column_vector& b) {
             	Eigen::MatrixXd kappaMat; kappaMat.setZero(numInteractions, N); 
             	dlibToEig(N, numInteractions, b, kappaMat);
             	return -computeFinalProb(N, num_states, problem, f_rate,
             	kappaMat, t_disc, initial, targets);}, 

             
             kappaVec,-0.99);

	return m;

}


void findProtocolSA_dlib_E(int np, Database* db, int initial, int target, bool useFile) {
	//determine a protocol to maximize prob of target state

	//get database info
	int num_states = db->getNumStates(); 

	//fill in the lumpMap for the DB
	printf("Getting indices of lumped permutation states...\n");
	lumpPerms(db);
	//get total number of lumped states
	int lumped = 0;
	for (int i = 0; i < num_states; i++) {
		int lump_name = db->lumpMap[i];
		if (lump_name > lumped) {
			lumped = lump_name;
		}
	}
	lumped++; //there are 1 more states than max index
	printf("There are %d lumped states.\n\n", lumped);

	//set up particle identity
	int* particleTypes = new int[np];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(np, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(np, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;


	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries
	Eigen::MatrixXd f_rate; f_rate.setZero(num_states,num_states);

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < num_states*num_states; i++) f_rate(i) = Tconst[i];
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(np, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//get the particle permutation
	readKappaFile(numInteractions, kappaVals);

	//set problem parameters
	double T = 5;                        //final time
	double rho = RANGE;                   //range parameter to morse potential

	//discretization parameters
	int N = 45;                            //number of intervals in termporal discretization
	double dt = double(T) / double(N);    //size of time intervals

	Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
	t_disc.setLinSpaced(N+1,0,T);
	//std::cout << t_disc << "\n";
	//abort();

	Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
	t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
	t_shift += t_disc;

	Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
	t_centers = t_shift.segment(0,N) / 2.0;

	//create discretized kappa vals
	Eigen::MatrixXd kappaDisc; kappaDisc.setZero(numInteractions, N);
	for (int i = 0; i < numInteractions; i++) {
		for (int j = 0; j < N; j++) {
			kappaDisc(i,j) = kappaVals[i];
		}
	}

	//set up dlib storage types
	column_vector kappaProt(N*numInteractions);
	eigToDlib(N, numInteractions, kappaDisc, kappaProt);

	//create an optInfo object for the problem
	OptInfo problem = OptInfo(np, num_states, numTypes, numInteractions, particleTypes, 
														initial, target, rho, db);

	//perform optimization - try to escape local minima
	bool same = false; double obj = 0; double tol = 1e-6; 
	int attempt = 0; int max_attempts = 1;
	while (!same) { //if two consecutive ptimization dont give same obj value, perturb
		double m = dlib_objective_call_fd(N, num_states, numInteractions, &problem, kappaProt, f_rate, t_disc, 
																	 t_centers, initial, targets);
		if (fabs(obj-m) < tol) { //escape
			same = true;
		}
		else { //perturb
			if (attempt > 0) {
				kappaProt /= 1.001;
				obj = m;
			}
		}
		attempt++;
		if (attempt >= max_attempts) {
			break;
		}
	}

	//go from column vec to eigen vec
	dlibToEig(N, numInteractions, kappaProt, kappaDisc);

	//print the found protocol to terminal
	std::cout << kappaProt << "\n";
	
	//output the protocol and probbility info to a file for matlab figure generation

	//create a matrix to store state probs as fn of time
	Eigen::MatrixXd probs; probs.setZero(N, num_states);
	computeProbVec(N, num_states, &problem, f_rate, kappaDisc, t_disc, 
															initial, targets, probs);

	//collapse the full state space onto lumped states
	Eigen::MatrixXd lumpProbs; lumpProbs.setZero(N, lumped);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < num_states; j++) {
			lumpProbs(i,db->lumpMap[j]) += probs(i,j);
		}
	}

	//initial probability vec
	Eigen::VectorXd initProb; initProb.setZero(lumped);
	initProb(db->lumpMap[initial]) = 1.0;

	//get the target probability
	int targetLump = db->lumpMap[target];
	double p = lumpProbs(N-1, targetLump);

	printf("Using the protocol \n\n");
	std::cout << kappaDisc << "\n\n";
	printf("We get a final probability of %f\n", p);

	std::cout << lumpProbs << "\n";

	std::ofstream ofile;
	ofile.open("non_eq/protocolTest.txt");
	outputProtocol(N, T, targetLump, initProb, lumpProbs, kappaDisc, ofile);
	ofile.close();


	//free memory
	delete []Tconst; delete []particleTypes; delete []kappaVals;


}

void findProtocolSA_adjoint(int np, Database* db, int initial, int target, bool useFile) {
	//determine a protocol to maximize prob of target state

	//get database info
	int num_states = db->getNumStates(); 

	//fill in the lumpMap for the DB
	printf("Getting indices of lumped permutation states...\n");
	lumpPerms(db);
	//get total number of lumped states
	int lumped = 0;
	for (int i = 0; i < num_states; i++) {
		int lump_name = db->lumpMap[i];
		if (lump_name > lumped) {
			lumped = lump_name;
		}
	}
	lumped++; //there are 1 more states than max index
	printf("There are %d lumped states.\n\n", lumped);

	//set up particle identity
	int* particleTypes = new int[np];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(np, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(np, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;


	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries
	Eigen::MatrixXd f_rate; f_rate.setZero(num_states,num_states);

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < num_states*num_states; i++) f_rate(i) = Tconst[i];
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(np, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//get the particle permutation
	readKappaFile(numInteractions, kappaVals);

	//set problem parameters
	double T = 5;                        //final time
	double rho = RANGE;                   //range parameter to morse potential

	//discretization parameters
	int N = 30;                            //number of intervals in termporal discretization
	double dt = double(T) / double(N);    //size of time intervals

	Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
	t_disc.setLinSpaced(N+1,0,T);
	//std::cout << t_disc << "\n";
	//abort();

	Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
	t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
	t_shift += t_disc;

	Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
	t_centers = t_shift.segment(0,N) / 2.0;

	//create discretized kappa vals
	Eigen::MatrixXd kappaDisc; kappaDisc.setZero(numInteractions, N);
	for (int i = 0; i < numInteractions; i++) {
		for (int j = 0; j < N; j++) {
			kappaDisc(i,j) = kappaVals[i];
		}
	}

	//create an optInfo object for the problem
	OptInfo problem = OptInfo(np, num_states, numTypes, numInteractions, particleTypes, 
														initial, target, rho, db);

	//do the optimization
	protocolOptAdjoint_energy(N, kappaDisc, num_states, &problem, f_rate, t_disc, t_centers, 
										 initial, targets);

	

	//print the found protocol to terminal
	std::cout << kappaDisc << "\n";
	
	//output the protocol and probbility info to a file for matlab figure generation

	//create a matrix to store state probs as fn of time
	Eigen::MatrixXd probs; probs.setZero(N, num_states);
	computeProbVec(N, num_states, &problem, f_rate, kappaDisc, t_disc, 
															initial, targets, probs);

	//collapse the full state space onto lumped states
	Eigen::MatrixXd lumpProbs; lumpProbs.setZero(N, lumped);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < num_states; j++) {
			lumpProbs(i,db->lumpMap[j]) += probs(i,j);
		}
	}

	//initial probability vec
	Eigen::VectorXd initProb; initProb.setZero(lumped);
	initProb(db->lumpMap[initial]) = 1.0;

	//get the target probability
	int targetLump = db->lumpMap[target];
	double p = lumpProbs(N-1, targetLump);

	printf("Using the protocol \n\n");
	std::cout << kappaDisc << "\n\n";
	printf("We get a final probability of %f\n", p);

	std::cout << lumpProbs << "\n";

	std::ofstream ofile;
	ofile.open("non_eq/protocolTest.txt");
	outputProtocol(N, T, targetLump, initProb, lumpProbs, kappaDisc, ofile);
	ofile.close();


	//free memory
	delete []Tconst; delete []particleTypes; delete []kappaVals;


}


void setProtocol(int N, Eigen::MatrixXd& protocol, int which) {
	//set example protocols

	for (int i = 0; i < N; i++) {

		double t = double(i) / double(N);

		if (which == 0) {
			protocol(0,i) = 10000;
			protocol(1,i) = 5 + t * 80.0;
			protocol(2,i) = 0.1 + t * 0.0;
		}

		if (which == 1) {
			protocol(0,i) = 10000;
			protocol(1,i) = 20;
			protocol(2,i) = 0.1 + t * 0.0;
		}
	}
}


void testProtocol_E(int np, Database* db, int initial, int target, bool useFile) {
	//test a supplied protocol
	//outputs the probability of each state as a function of time

	//get database info
	int num_states = db->getNumStates(); 

	//fill in the lumpMap for the DB
	printf("Getting indices of lumped permutation states...\n");
	lumpPerms(db);
	//get total number of lumped states
	int lumped = 0;
	for (int i = 0; i < num_states; i++) {
		int lump_name = db->lumpMap[i];
		if (lump_name > lumped) {
			lumped = lump_name;
		}
	}
	lumped++; //there are 1 more states than max index
	printf("There are %d lumped states.\n\n", lumped);


	//set up particle identity
	int* particleTypes = new int[np];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(np, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(np, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;


	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries
	Eigen::MatrixXd f_rate; f_rate.setZero(num_states,num_states);

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < num_states*num_states; i++) f_rate(i) = 0.01*Tconst[i];
	printf("All the ground states are:\n");
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(np, num_states, target, db, targets);
	printf("All the target states are:\n");
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//get the particle permutation
	readKappaFile(numInteractions, kappaVals);

	//set problem parameters
	double T = 200;                        //final time
	double rho = RANGE;                   //range parameter to morse potential

	//discretization parameters
	int N = 45;                            //number of intervals in termporal discretization
	double dt = double(T) / double(N);    //size of time intervals

	Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
	t_disc.setLinSpaced(N+1,0,T);

	Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
	t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
	t_shift += t_disc;

	Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
	t_centers = t_shift.segment(0,N) / 2.0;


	//************************************************************************//
	//create discretized kappa vals - input test values here 
	Eigen::MatrixXd kappaDisc; kappaDisc.setZero(numInteractions, N);
	for (int i = 0; i < numInteractions; i++) {
		for (int j = 0; j < N; j++) {
			kappaDisc(i,j) = kappaVals[i];
		}
	}
	for (int j = 0; j < N; j++) {
		if (j < 30) {
			kappaDisc(0,j) = 0.1;
			//kappaDisc(1,j) = 0.1;
			//kappaDisc(2,j) = 0.1;
		}
		else {
			kappaDisc(0,j) = 10000;
			//kappaDisc(1,j) = 10000;
			//kappaDisc(2,j) = 10000;
		}
	}
	//********************************************************************//

	//make a special protocol
	setProtocol(N, kappaDisc, 1);

	//create an optInfo object for the problem
	OptInfo problem = OptInfo(np, num_states, numTypes, numInteractions, particleTypes, 
														initial, target, rho, db);


	//create a matrix to store state probs as fn of time
	Eigen::MatrixXd probs; probs.setZero(N, num_states);
	computeProbVec(N, num_states, &problem, f_rate, kappaDisc, t_disc, 
															initial, targets, probs);

	//collapse the full state space onto lumped states
	Eigen::MatrixXd lumpProbs; lumpProbs.setZero(N, lumped);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < num_states; j++) {
			lumpProbs(i,db->lumpMap[j]) += probs(i,j);
		}
	}

	//initial probability vec
	Eigen::VectorXd initProb; initProb.setZero(lumped);
	initProb(db->lumpMap[initial]) = 1.0;

	//get the target probability
	int targetLump = db->lumpMap[target];
	double p = lumpProbs(N-1, targetLump);

	printf("Using the protocol \n\n");
	std::cout << kappaDisc << "\n\n";
	printf("We get a final probability of %f\n", p);

	std::cout << lumpProbs << "\n";

	std::ofstream ofile;
	ofile.open("non_eq/protocolTest.txt");
	outputProtocol(N, T, targetLump, initProb, lumpProbs, kappaDisc, ofile);
	ofile.close();

	delete []particleTypes; delete []kappaVals; delete []Tconst;




}

void outputProtocol(int N, double T, int targetLump, const Eigen::VectorXd& initProb, 
										const Eigen::MatrixXd& probs, const Eigen::MatrixXd& prot, 
										std::ostream& ofile) {
	//output the protocol data as well as the probability of each state
	//as a function of time

	ofile << N << "\n";
	ofile << T << "\n";
	ofile << targetLump << "\n";
	ofile << prot << "\n";
	ofile << initProb.transpose() << "\n";
	ofile << probs << "\n";
}






/* testing */

void testEqDeriv(int N, Database* db) {
	//test if the analytic derivative of eq distr. is correct.
	//test result - correct!

	//get database info
	int num_states = db->getNumStates(); 

	bool useFile = true;
	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//get the particle permutation
	readKappaFile(numInteractions, kappaVals);
	Eigen::VectorXd kappaVec; kappaVec.setZero(numInteractions);
	for (int i = 0; i < numInteractions; i++) kappaVec(i) = kappaVals[i];

	//set other parameters to call eq gradient
	int bond = 1;
	double* eq = new double[num_states]; double* eqDeriv = new double[num_states];
	for (int i = 0; i < num_states; i++) eq[i] = 0;
	for (int i = 0; i < num_states; i++) eqDeriv[i] = 0;


	//eval the eq prob at the given kappa
	std::map<std::pair<int,int>,double> kappaMap;    //map for interaction type to kappa
	makeKappaMap(numTypes, kappaVec, kappaMap);
	reweight(N, num_states, db, particleTypes, eq, kappaMap);

	//use the eq to compute the analytical gradient
	eqDerivative(N, num_states, numTypes, numInteractions, db, bond, eq, eqDeriv, 
		particleTypes, kappaVec);

	//eval at kappa+h for fd
	double h = 1e-6; kappaVec[bond] += h;
	double* eq2 = new double [num_states]; for (int i = 0; i < num_states; i++) eq2[i] = 0;
	double* fd = new double [num_states]; for (int i = 0; i < num_states; i++) fd[i] = 0;
	kappaMap.clear(); makeKappaMap(numTypes, kappaVec, kappaMap);
	reweight(N, num_states, db, particleTypes, eq2, kappaMap);

	for (int i = 0; i < num_states; i++) {
		fd[i] = (eq2[i] - eq[i]) / h;
	}

	for (int i = 0; i < num_states; i++) {
		printf("state %d, analytical %f, fd %f\n", i, eqDeriv[i], fd[i]);
	}

	delete []eq; delete []eq2; delete []eqDeriv; delete []fd;


}


void testProbGradient(int np, Database* db) {
	//test if the gradient calculation is the same as fd of the probability fn

	int initial = 1; int target = 10; bool useFile = true;

	//get database info
	int num_states = db->getNumStates(); 

	//fill in the lumpMap for the DB
	printf("Getting indices of lumped permutation states...\n");
	lumpPerms(db);
	//get total number of lumped states
	int lumped = 0;
	for (int i = 0; i < num_states; i++) {
		int lump_name = db->lumpMap[i];
		if (lump_name > lumped) {
			lumped = lump_name;
		}
	}
	lumped++; //there are 1 more states than max index
	printf("There are %d lumped states.\n\n", lumped);

	//set up particle identity
	int* particleTypes = new int[np];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = readDesignFile(np, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = setTypes(np, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;


	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries
	Eigen::MatrixXd f_rate; f_rate.setZero(num_states,num_states);

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < num_states*num_states; i++) f_rate(i) = 0.2*Tconst[i];
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	findIsomorphic(np, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//get the particle permutation
	readKappaFile(numInteractions, kappaVals);

	//set problem parameters
	double T = 1;                        //final time
	double rho = RANGE;                   //range parameter to morse potential

	//discretization parameters
	int N = 20;                            //number of intervals in termporal discretization
	double dt = double(T) / double(N);    //size of time intervals

	Eigen::VectorXd t_disc; t_disc.setZero(N+1);    //temporal discretization
	t_disc.setLinSpaced(N+1,0,T);

	Eigen::VectorXd t_shift; t_shift.setZero(N+1);  //temporal shift vector
	t_shift.segment(0,N) = t_disc.segment(1,N); t_shift(N) = t_disc(0);
	t_shift += t_disc;

	Eigen::VectorXd t_centers; t_centers.setZero(N); //midpoints of temporal intervals
	t_centers = t_shift.segment(0,N) / 2.0;

	//create discretized kappa vals
	Eigen::MatrixXd kappaDisc; kappaDisc.setZero(numInteractions, N);
	for (int i = 0; i < numInteractions; i++) {
		for (int j = 0; j < N; j++) {
			kappaDisc(i,j) = kappaVals[i];
		}
	}

	//create an optInfo object for the problem
	OptInfo problem = OptInfo(np, num_states, numTypes, numInteractions, particleTypes, 
														initial, target, rho, db);

	//set up storage for the gradients
	Eigen::VectorXd R; R.setZero(N*numInteractions);
	Eigen::VectorXd Rfd; Rfd.setZero(N*numInteractions);
	Eigen::VectorXd Ra; Ra.setZero(N*numInteractions);


	//compute the gradients
	evalGradientAdjoint(N, kappaDisc, num_states, &problem, f_rate, t_disc, 
							 t_centers, initial, targets, Ra);

	printf("Adjoint method finished\n");

	evalGradientFD(N, kappaDisc, num_states, &problem, f_rate, t_disc, 
							 t_centers, initial, targets, Rfd);

	printf("FD Method finished\n");

	//evalGradient(N, kappaDisc, num_states, &problem, f_rate, t_disc, 
	//						 t_centers, initial, targets, R);

	printf("Primal method finished\n");



	for (int i = 0; i < N*numInteractions; i++) {
		printf("interval %d, analytical %f, fd %f, ratio %f\n", i % N, R(i), Rfd(i), R(i)/Rfd(i));
	}

	for (int i = 0; i < N*numInteractions; i++) {
		printf("interval %d, adjoint %f, fd %f, ratio %f\n", i % N, Ra(i), Rfd(i), Ra(i)/Rfd(i));
	}


	delete []particleTypes; delete []Tconst;
	delete []kappaVals;



}




















}