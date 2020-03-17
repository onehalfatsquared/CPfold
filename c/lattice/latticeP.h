#pragma once
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <utility>
#include <random>
#include <chrono>


struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const
    {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

class RandomNo{
	std::mt19937 generator; 
	std::uniform_real_distribution<double> uDist;
	std::normal_distribution<double> gDist;

	public:
		RandomNo(unsigned int seed = std::random_device{}())
        : generator{seed},uDist{0,1},gDist{0,1} {}

    double getU() {
    	return uDist(generator);
    }
    double getG() {
    	return gDist(generator);
    }
};

namespace lattice {


/* Particle struct stores lattice coordinates and type */
struct Particle{
	//null constructor
	Particle() {x = y = type = 0;}
		
	//declare variables
	int type;
	int x; int y;

	//non-null constructor
	Particle(int x_, int y_, int type_) :  x(x_), y(y_), type(type_) {}
	Particle(int x_, int y_) :  x(x_), y(y_), type(0) {}
};

typedef std::unordered_map<std::pair<int,int>,Particle*,pair_hash> particleMap; 


//functions to set up a chain
void getTypes(int N, int* types, bool useFile);
void initChain(int N, Particle* chain, particleMap& cMap, 
							 bool useFile);





//functions to determine moves
void checkRotation(int x, int y, int particle, Particle* chain, std::vector<std::pair<int,int>>& moves,
									 particleMap& cMap);

void getRotations(int particle, int neighbor, Particle* chain, std::vector<std::pair<int,int>>& moves,
						  particleMap& cMap);

void getCorners(int particle, Particle* chain, std::vector<std::pair<int,int>>& moves,
						  particleMap& cMap);

void getMoves(int N, int particle, Particle* chain, std::vector<std::pair<int,int>>& moves,
						  particleMap& cMap);


//energy functions
void getBonds(int N, Particle* chain, std::vector<std::pair<int,int>> bonds);

//mcmc functions
int randomInteger(int N, RandomNo* rngee);
void acceptMove(int particle, int x_old, int y_old, Particle* chain, 
								particleMap& cMap);

void rejectMove(int particle, int x_old, int y_old, Particle* chain);
void takeStep(int N, Particle* chain, particleMap& cMap,
							RandomNo* rngee, double& energy);



}