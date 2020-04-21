#pragma once
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <utility>
#include <random>
#include <chrono>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include "pair.h"
#include "design.h"
#include "nauty.h"
#include "../defines.h"

namespace ga {



class Person {
public:

	//con(de)structors
	Person();
	~Person();

	Person(int num_interactions, double* kappaVals);

	void copy(const Person& old);
	Person(const Person& old) {
		copy(old);
	}

	void destroy();

	Person& operator=(const Person& old) {
		if (this != &old) {
			destroy();
			copy(old);
		}
		return *this;
	}

	//parameters
	int num_interactions;
	double Rate;
	double Eq;
	double fitness;
	double* kappaVals;
	std::map<std::pair<int,int>,double> kappa;

	//functions 
	void setKappa(int num_interactions, double* kappaVals);
	Person mate(Person partner);
	void evalStats(int N, bd::Database* db, int initial, std::vector<int> targets, 
								 int* particleTypes, double* eq, double* Tconst, double* T, double* m);
	void evalFitness(double eq, double rate);

private:








};





void perform_evolution(int N, bd::Database* db, int initial, int target, bool useFile);























}