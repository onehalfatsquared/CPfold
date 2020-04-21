#include "genetics.h"

#define PI 3.1415926

namespace ga {


Person::Person() {
	Rate = Eq = 0;
	num_interactions = 0;
	kappaVals = NULL;
}

Person::~Person() {
	destroy();
}

void Person::copy(const Person& old) {
	Rate = old.Rate;
	Eq = old.Eq;
	num_interactions = old.num_interactions;

	kappaVals = new double[num_interactions];
	for (int i = 0; i < num_interactions; i++) {
		kappaVals[i] = old.kappaVals[i];
	}

	kappa = old.kappa;
}

void Person::destroy() {
	delete []kappaVals;
}

Person::Person(int num_interactions_, double* kV) {
	Rate = Eq = 0;
	num_interactions = num_interactions_;
	kappaVals = new double[num_interactions];
	for (int i = 0; i < num_interactions; i++) {
		kappaVals[i] = kV[i];
	}
	bd::makeKappaMap(2, kappaVals, kappa);
}

void Person::setKappa(int num_interactions_, double* kV) {
	//construct the array and set the kappa values 
	num_interactions = num_interactions_;
	kappaVals = new double[num_interactions];
	for (int i = 0; i < num_interactions; i++) {
		kappaVals[i] = kV[i];
	}
	bd::makeKappaMap(2, kappaVals, kappa);
}

Person Person::mate(Person partner) {
	//create offspring from two parents

	//create new trait array
	double* kV = new double[num_interactions];

	//take each trait w/ 50/50 chance from each parent
	for (int i = 0; i < num_interactions; i++) {
		double p = double(rand()) / RAND_MAX;
		
		if (p < 0.45) { //parent 1
			kV[i] = kappaVals[i];
		}
		else if (p > 0.45 && p < 0.9) { //parent 2
			kV[i] = partner.kappaVals[i];
		}
		else { //mutate
			kV[i] = tan(PI/2.0*rand()/RAND_MAX);
		}

		if (kV[i] > 100 && kV[i] < 1e4) {
			double p2 = double(rand()) / RAND_MAX;
			if (p2 < 0.45) {
				kV[i] *= 10;
			}
		}

		if (kV[i] > 1e3) {
			kV[i] = 1e5;
		}

	}

	//create the offspring
	Person kid = Person(num_interactions, kV);

	//free memory
	delete []kV;

	//birth the child
	return kid;
}

Person createNew(int num_interactions) {
	//create a new person w/o mating

	//create new trait array
	double* kV = new double[num_interactions];
	for (int j = 0; j < num_interactions; j++) {
		kV[j] = tan(PI/2.0*rand()/RAND_MAX);
	}

	//create a person and return it
	Person immigrant = Person(num_interactions, kV);
	return immigrant;
}

void Person::evalStats(int N, bd::Database* db, int initial, std::vector<int> targets, 
											 int* particleTypes, double* eq, double* Tconst, double* T, double* m) {
	//compute the eq and rate for this person

	int num_states = db->getNumStates();

	//do rewieght
	bd::reweight(N, num_states, db, particleTypes, eq, kappa);
	//get eq prob
	double eqProb = bd::getEqProb(initial, kappaVals, db, particleTypes, targets);
	//copy Tconst into T
	std::copy(Tconst, Tconst+num_states*num_states, T);
	//fill in transposed entries such that T satisfies detailed balance
	bd::satisfyDB(T, num_states, db, eq);
	//fill in diagonal with negative sum of row entries
	bd::fillDiag(T, num_states);
	//get the transition rate
	bd::computeMFPTs(num_states, T, targets, m);
	double rate = 1/m[initial];

	Rate = rate; Eq = eqProb;
}

void Person::evalFitness(double eqMax, double rateMax) {
	//evulate the fitness as weighted sum of eq and rate

	fitness = Eq / eqMax + Rate / rateMax;
}


bool is_dominated(int pop_size, double x, double y, double* xAll, double* yAll) {
	//determine of the points x and y are dominated by any others in the set
	double tol = -1e-6;

	for (int i = 0; i < pop_size; i++) {
		double d1 = x - xAll[i]; double d2 = y - yAll[i];
		bool s1 = (d1 > tol); bool s2 = (d2 > tol);
		//printf("Diffs %f and %f, bools %d and %d\n", d1, d2, s1, s2);
		if (!s1 && !s2) {
			return true;
		}
	}

	return false;
}

void non_dominated_set(int pop_size, double* xAll, double* yAll, std::vector<int>& nonDom) {
	//find the elements of xAll and yAll that are non-dominated by anything else in the set

	for (int i = 0; i < pop_size; i++) {
		//printf("DOing person %d\n\n\n", i);
		double x = xAll[i]; double y = yAll[i];
		bool D = is_dominated(pop_size, x, y, xAll, yAll);
		if (!D) {
			nonDom.push_back(i);
		}
	}
}









void perform_evolution(int N, bd::Database* db, int initial, int target, bool useFile) {
	//run the genetic algorithm on the colloid database to find the pareto front
	//for the target state. Uses fixed particle types from file. 

	//parameters to the genetic algorithm
	int generations = 400;
	int pop_size    = 1500;
	double elite_p  = 0.1;
	double mates_p  = 0.5;

	//if we have prior estimates of the max rate and eqProb, set here.
	//otherwise, this will update, adaptively. 
	double rateMax = 0.5; double eqMax = 0.8;

	//get database info
	int num_states = db->getNumStates(); 

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = bd::readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		int IC = 1; 
		numTypes = bd::setTypes(N, particleTypes, IC);
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//declare rate matrix, probability transition matrix, equilibrium measure
	double* T = new double[num_states*num_states]; //rate matrix
	double* Tconst = new double[num_states*num_states]; //rate matrix - only forward entries
	double* eq = new double[num_states];           //equilibrium measure
	double* m = new double[num_states];            //mfpts 

	//init the rate matrix with zeros
	for (int i = 0; i < num_states*num_states; i++) {
		Tconst[i] = 0;
	}

	//get bonds->bonds+1 entries from mfpt estimates
	std::vector<int> ground; //vector to hold all ground states
	bd::createTransitionMatrix(Tconst, num_states, db, ground);
	for (int i = 0; i < ground.size(); i++) {
		std::cout << ground[i] << "\n";
	}

	//find all target states consistent with input target
	std::vector<int> targets; 
	bd::findIsomorphic(N, num_states, target, db, targets);
	for (int i = 0; i < targets.size(); i++) {
		std::cout << targets[i] << "\n";
	}

	//construct the initial population - tan(uniformly random)
	std::vector<Person> population;
	for (int i = 0; i < pop_size; i++) {
		for (int j = 0; j < numInteractions; j++) {
			kappaVals[j] = tan(PI/2.0*rand()/RAND_MAX);
		}
		Person p = Person(numInteractions, kappaVals);
		p.evalStats(N, db, initial, targets, particleTypes, 
							  eq, Tconst, T, m);
		p.evalFitness(eqMax, rateMax);
		population.push_back(p);
	}

	//declare outfile
	std::ofstream ofile;
	ofile.open("paretoGA.txt");

	//init storage for storing eq and rate
	double* popEq = new double[pop_size];
	double* popRate = new double[pop_size];
	double* popFitness = new double[pop_size];

	//loop over generations
	for (int gen = 0; gen < generations; gen++) {
		//get the eq, rate, and fitness of each person

		printf("Spawning generation %d of %d\n", gen+1,generations);
		for (int i = 0; i < pop_size; i++) {
			popEq[i] = population[i].Eq; 
			popRate[i] = population[i].Rate; 
			popFitness[i] = population[i].fitness;

			//printf("Person %d, eq %f, rate %f\n", i, population[i].Eq, population[i].Rate);
			//printf("%f, %f, %f\n", population[i].kappaVals[0],population[i].kappaVals[1],population[i].kappaVals[2]);
		}

		//get the current max eq and rate to update the scalings
		double eM = *std::max_element(popEq, popEq + pop_size); 
		double rM = *std::max_element(popRate, popRate + pop_size); 
		if (eM > eqMax) {
			eqMax = eM;
		}
		if (rM > rateMax) {
			rateMax = rM;
		}

		//next, we sort the population by fitness, high to low
		std::vector<int> p(pop_size);
		std::iota(p.begin(), p.end(), 0);
		std::sort(p.begin(), p.end(), [&](int i1, int i2) { return popFitness[i1] > popFitness[i2]; });
		
		//determine the non-dominated points
		std::vector<int> nonDom;
		non_dominated_set(pop_size, popEq, popRate, nonDom);
		
		//create the new generation, perform elitism step
		std::vector<Person> new_generation;
		int elites = nonDom.size(); 
		for (int i = 0; i < elites; i++) {
			new_generation.push_back(population[nonDom[i]]);
		}

		//fill rest by mating the top percent of the present generation
		int rest = pop_size - elites;
		int top = mates_p * pop_size;
		for (int i = 0; i < rest; i++) {
			int r1 = rand() % top;
			int r2 = rand() % top;
			Person p1 = population[r1];
			Person p2 = population[r2];
			Person kid = p1.mate(p2);
			kid.evalStats(N, db, initial, targets, particleTypes, 
							  	eq, Tconst, T, m);
			kid.evalFitness(eqMax, rateMax);
			new_generation.push_back(kid);
		}

		population = new_generation;

	}

	//output the final results
	for (int i = 0; i < pop_size; i++) {
		ofile << population[i].Eq << ' ' << population[i].Rate << ' ';
		ofile << population[i].kappaVals[0] << ' ';
		ofile << population[i].kappaVals[1] << ' ';
		ofile << population[i].kappaVals[2] << ' ' << "\n";;
	}

	ofile.close();





	//free memory
	delete []particleTypes; delete []kappaVals; delete []T; delete []Tconst;
	delete []eq; delete []m;
	delete []popRate; delete []popEq; delete []popFitness;


}















}