#include "genetics.h"
#include "bDynamics.h"

#define PI 3.1415926

namespace ga {


void readTargetFile(int N, int* M_target) {
	//read in bonds present in the target state

	//file location
	std::string filename = "input/design/target.txt";

	//open the file
	std::ifstream in_str(filename);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		return;
	}

	//store the entries
	int p1; int p2;
	while (in_str >> p1) {
		//read line by line, add entry into AM
		in_str >> p2;

		M_target[bd::toIndex(p1,p2,N)] = 1;
		M_target[bd::toIndex(p2,p1,N)] = 1;
	}

	//fill in super and sub diagonal
	for (int i = 0; i < N-1; i++) {
		M_target[bd::toIndex(i,i+1,N)] = 1;
		M_target[bd::toIndex(i+1,i,N)] = 1;
	}
}

void doSampling(int N, double Tf, int* types, std::map<std::pair<int,int>,double> kappa,
							  double& eq, double& time, int* M_target) {
	//sample a trajectory and compute the fraction of time spent in target as well
	//as the minimum time the target was hit

	//init a chain
	double* X = new double[DIMENSION*N];
	int* M = new int[N*N]; for (int i = 0; i < N*N; i++) M[i] = 0;
	bd::setupChain(X,N);

	bool firstHit = true;
	double firstHitTime = 0;
	double eqTime = 0;

	//set parameters
	double t = 0; double DT = 0.01;
	double beta = 1.0;  double rho = 10;
	double method = 1; double pot = 0;
	double cut = 1.1;

	int* P = new int[N*N];
	double* E = new double[N*N];
	bd::fillP(N, types, P, E, kappa);

	while (t < Tf) {
		bd::solveSDE(X, N, DT, rho, beta, E, P, method, pot);
		t += DT;

		bd::getAdjCut(X,N,M,cut);
		/*
		bd::printAM(N,M);
		std::cout << "\n";
		bd::printAM(N,M_target);
		abort();
		*/
		int same = bd::checkSame(M, M_target, N);

		if (same) {
			eqTime += DT;
			//printf("HIT\n");
			if (firstHit) {
				firstHitTime = t;
				firstHit = false;
			}
		}
	}

	eq = eqTime / Tf; time = firstHitTime;
}

void Person::evalStats(double Tf, int samples, int* M_target) {
	//evaluate a rate and eq prob surrogate by sampling

	//estimate objectives from given number of sample trajectories
	double eqEstimate = 0; double rateEstimate = 0;

	//get the samples
	for (int i = 0; i < samples; i++) {
		//get sample
		double eqSample = 0; double timeSample = 0;
		doSampling(N, Tf, types, kappa, eqSample, timeSample, M_target);

		//update eq time estimate
		eqEstimate += eqSample; 

		//update rate estimate. gets 0 if structure did not form
		if (timeSample > 0.001) {
			rateEstimate += 1.0 / timeSample;
		}
	}

	//weight by num samples
	eqEstimate /= double(samples);
	rateEstimate /= double(samples);

	//set the class variables
	Eq = eqEstimate;
	Rate = rateEstimate;

} 

void Person::applyBound(double lower, double upper) {
	//apply bounds to kappa values

	for (int i = 0; i < num_interactions; i++) {
		if (kappaVals[i] > upper) {
			kappaVals[i] = upper;
		}
		else if (kappaVals[i] < lower) {
			kappaVals[i] = lower;
		}
	}

	bd::makeKappaMap(numTypes, kappaVals, kappa);
}






void perform_evolution_sampling(int N, bool useFile) {
	//run the genetic algorithm on the colloid database to find the pareto front
	//for the target state. Uses fixed particle types from file. 

	//get rid of the eigen parallelism
	Eigen::setNbThreads(0);

	//parameters to the genetic algorithm
	int generations = 100;
	int pop_size    = 70;
	double elite_p  = 0.1;
	double mates_p  = 0.2;
	bool printAll   = false;         //set true to make movie of output

	//sampling parameters
	double Tf = 5.0;
	int samples = 40;

	//if we have prior estimates of the max rate and eqProb, set here.
	//otherwise, this will update, adaptively. 
	double rateMax = 0.1; double eqMax = 0.1;

	//set up particle identity
	int* particleTypes = new int[N];
	int numTypes;
	if (useFile) { //use the fle to set identities
		numTypes = bd::readDesignFile(N, particleTypes);
	}
	else { //uses the function to set identities
		//int IC = 1; 
		//numTypes = bd::setTypes(N, particleTypes, IC);
		numTypes = 3;
	}
	int numInteractions = numTypes*(numTypes+1)/2;

	//set up sticky parameter values
	double* kappaVals = new double[numInteractions];

	//set up target
	int* M_target = new int[N*N]; for (int i = 0; i < N*N; i++) M_target[i] = 0;
	readTargetFile(N, M_target);

	//construct the initial population
	printf("Generating the initial population\n");
	Person* pop_array = new Person[pop_size];
	std::vector<Person> population;
	#pragma omp parallel 
	{
	//declare all arrays we need to do calculations
	RandomNo* rngee = new RandomNo();              //random number generator

	//loop over population
	#pragma omp for schedule(dynamic)
	for (int i = 0; i < pop_size; i++) {
		//draw parameters from a distribution
		sampleParameters(N, numInteractions, kappaVals, particleTypes, numTypes, useFile, rngee);
		
		//create a person, evaluate their stats
		Person p = Person(N, numInteractions, numTypes, particleTypes, kappaVals);
		p.applyBound(0.1, 1000);
		p.evalStats(Tf, samples, M_target);
		p.evalFitness(eqMax, rateMax);
		pop_array[i] = p;
		//printf("e %f, r %f, f %f\n", pop_array[i].Eq, pop_array[i].Rate, pop_array[i].fitness);
		printf("Finsihing sample %d on thread %d\n", i, omp_get_thread_num());
	}
	//free memory
	delete rngee;
	//end parallel region
	}

	//move from array to vector
	for (int i = 0; i < pop_size; i++) {
		population.push_back(pop_array[i]);
	}
	delete []pop_array;

	//declare outfile
	std::ofstream ofile;
	ofile.open("paretoGAsampling.txt");

	//print the initial population
	if (printAll)
		printPopulation(population, pop_size, ofile);

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

		//printf("Gen %d, em %f, rm %f\n", gen, eqMax, rateMax);

		//next, we sort the population by fitness, high to low
		std::vector<int> p(pop_size);
		std::iota(p.begin(), p.end(), 0);
		//std::sort(p.begin(), p.end(), [&](int i1, int i2) { return popFitness[i1] > popFitness[i2]; });
		
		double f = population[p[0]].fitness;
		std::cout << "max f " << f << "\n";

		//determine the non-dominated points
		std::vector<int> nonDom;
		non_dominated_set(pop_size, popEq, popRate, nonDom);
		
		//create the new generation, perform elitism step
		std::vector<Person> new_generation;
		int elites = nonDom.size(); 
		printf("Gen %d, num elites %d\n", gen, elites);
		for (int i = 0; i < elites; i++) {
			new_generation.push_back(population[nonDom[i]]);
		}

		//fill rest by mating the top percent of the present generation
		int rest = pop_size - elites;
		int top = mates_p * pop_size;
		Person* pop_array = new Person[rest];
		#pragma omp parallel 
		{
		//init the arrays
		RandomNo* rngee = new RandomNo(); 
		//loop over population
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < rest; i++) {
			int r1 = p[floor(rngee->getU()*top)];
			int r2 = p[floor(rngee->getU()*top)];
			Person p1 = population[r1];
			Person p2 = population[r2];
			Person kid = p1.mate(p2, useFile, rngee);
			kid.applyBound(0.1, 1000);
			kid.evalStats(Tf, samples, M_target);
			kid.evalFitness(eqMax, rateMax);
			pop_array[i] = kid;
		}
		//end parallel region / free memory
		delete rngee;
		}

		//move from array to vector
		for (int i = 0; i < rest; i++) {
			new_generation.push_back(pop_array[i]);
		}
		delete []pop_array;

		//set the population equal to the newly generated one
		population = new_generation;

		//print if required
		if (printAll) {
			printPopulation(population, pop_size, ofile);
		}


	}

	//output the final results
	if (!printAll)
		printPopulation(population, pop_size, ofile);

	ofile.close();

	//print the particle types
	printTypes(population, pop_size, N);





	//free memory
	delete []particleTypes; delete []kappaVals; delete []M_target;
	delete []popRate; delete []popEq; delete []popFitness;


}








































}