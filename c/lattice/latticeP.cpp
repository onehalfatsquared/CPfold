#include "latticeP.h"



namespace lattice {


/******************************************************************************/
/***************** Setup and General Functions *******************************/
/******************************************************************************/


void getTypes(int N, int* types, bool useFile) {
	//fill in the types of particles

	if (useFile) {
		//read it from the file

	}
	else{
		//set to all zeros
		for (int i = 0; i < N; i++) {
			types[i] = 0;
		}
	}
	
}

void initChain(int N, Particle* chain, particleMap& cMap,
							 bool useFile) {
	//initialize a linear chain, with some type distribution

	//get the types of each particle
	int* types = new int[N];
	getTypes(N, types, useFile);

	//initialize the chain - linear on x axis
	for (int i = 0; i < N; i++) {
		chain[i] = Particle(i, 0, types[i]);
		cMap[std::make_pair(i, 0)] = &(chain[i]); 
	}

	//free type memory
	delete []types;

}

/******************************************************************************/
/***************** Monte Carlo Moves ******************************************/
/******************************************************************************/

void checkRotation(int x, int y, int particle, Particle* chain, std::vector<std::pair<int,int>>& moves,
									 particleMap& cMap) {
	//check if rotating to (x,y) violates any rules, add to moves if doesnt

	//first check if (x,y) is in the cMap
	particleMap::iterator occupied = cMap.find(std::make_pair(x,y));
	if (occupied != cMap.end()) { //position found in map
		return;
	}

	//check for 180 degree rotation
	int xP = chain[particle].x; int yP = chain[particle].y;
	if (abs(x-xP) == 2 || abs(y-yP) == 2) {
		return;
	}

	//if we reach here, this rotation can be added to moves
	moves.push_back(std::make_pair(x,y));
}

void getRotations(int particle, int neighbor, Particle* chain, std::vector<std::pair<int,int>>& moves,
						  particleMap& cMap) {
	//compute all valid rotations of an end particle

	//get the coordinates to rotate about
	int xN = chain[neighbor].x; int yN = chain[neighbor].y;

	//check all rotations
	checkRotation(xN-1, yN, particle, chain, moves, cMap);
	checkRotation(xN+1, yN, particle, chain, moves, cMap);
	checkRotation(xN, yN-1, particle, chain, moves, cMap);
	checkRotation(xN, yN+1, particle, chain, moves, cMap);

}

void getCorners(int particle, Particle* chain, std::vector<std::pair<int,int>>& moves,
						  particleMap& cMap) {
	//compute all valid corner moves

	//get coordinates of the neigboring particles
	int x = chain[particle].x;    int y = chain[particle].y;
	int xL = chain[particle-1].x; int yL = chain[particle-1].y; 
	int xR = chain[particle+1].x; int yR = chain[particle+1].y; 

	//check if this particle is at a corner - add the appropriate move
	if (abs(xL-xR) == 1 && abs(yL-yR) == 1) {
		int slope = (xL-xR) / (yL-yR);
		int pos = yL + slope * (x - xL);
		if (slope > 0) {
			if (y > pos) {
				moves.push_back(std::make_pair(x+1,y-1));
			}
			if (y < pos) {
				moves.push_back(std::make_pair(x-1,y+1));
			}
		}
		if (slope < 0) {
			if (y > pos) {
				moves.push_back(std::make_pair(x-1,y-1));
			}
			if (y < pos) {
				moves.push_back(std::make_pair(x+1,y+1));
			}
		}
	}

}

void getMoves(int N, int particle, Particle* chain, std::vector<std::pair<int,int>>& moves,
						  particleMap& cMap) {
	//get all moves for a given particle

	//check for end particles
	if (particle == 0) {
		getRotations(particle, particle+1, chain, moves, cMap);
	}
	else if (particle == N-1) {
		getRotations(particle, particle-1, chain, moves, cMap);
	}
	else { //this is an interior particle
		getCorners(particle, chain, moves, cMap);
	}

}

/******************************************************************************/
/***************** Energy Functions ******************************************/
/******************************************************************************/

void getBonds(int N, Particle* chain, std::vector<std::pair<int,int>> bonds) {
	//determine the non-trivial list of bonds

	for (int i = 0; i < N; i++) {
		int xi = chain[i].x; int yi = chain[i].y;

		for (int j = i+2; j < N; j++) {
			int xj = chain[j].x; int yj = chain[j].y;

			//check if particles are distance 1 apart. 1-norm distance
			int dist = abs(xj-xi) + abs(yj-yi);
			if (dist == 1) {
				bonds.push_back(std::make_pair(i,j));
			}
		}
	}

}


/******************************************************************************/
/***************** Monte Carlo Functions **************************************/
/******************************************************************************/

int randomInteger(int N, RandomNo* rngee) {
	//pick an integer uniformly from 0 to N-1

	double U = N * rngee->getU();
	return floor(U);

}

void acceptMove(int particle, int x_old, int y_old, Particle* chain, 
								particleMap& cMap) {
	//replace the old cMap entry with the new one
	cMap.erase(std::make_pair(x_old,y_old));

	int x_new = chain[particle].x; int y_new = chain[particle].y;
	cMap[std::make_pair(x_new, y_new)] = &(chain[particle]);

}

void rejectMove(int particle, int x_old, int y_old, Particle* chain) {
	//put the old coordinates back into the correct chain entry

	chain[particle].x = x_old; chain[particle].y = y_old;

}

void takeStep(int N, Particle* chain, particleMap& cMap,
							RandomNo* rngee, double& energy) {
	//perform an MCMC step - return the energy of the returned state

	//set the initial energy
	double e0 = energy;

	//first we pick a random particle and get its coordinates
	int particle = randomInteger(N, rngee);
	int x_old = chain[particle].x; int y_old = chain[particle].y;

	//next we generate the set of moves
	std::vector<std::pair<int,int>> moves;
	getMoves(N, particle, chain, moves, cMap);

	//pick a move to perform, or return if there are no choices
	int M = moves.size();
	if (M == 0) {
		return;
	}
	int move = randomInteger(M, rngee);

	//update the position of particle in chain to whats in move
	chain[particle].x = moves[move].first; chain[particle].y = moves[move].second; 

	//get energy from chain
	std::vector<std::pair<int,int>> bonds;
	getBonds(N, chain, bonds);
	double e1 = 2*bonds.size();

	//get acc probability - do accept/reject step
	double a = std::min(1.0, exp(-(e1-e0)));
	double U = rngee->getU();
	if (U <= a) {
		acceptMove(particle, x_old, y_old, chain, cMap);
		energy = e1;
	}
	else {
		rejectMove(particle, x_old, y_old, chain);
	}


}




}