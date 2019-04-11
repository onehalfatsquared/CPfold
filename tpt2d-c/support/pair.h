#pragma once
#include <cstdlib>

namespace bd{


/* pair structure that will store an index and a value
	 to be used to sparsely store transition probabilities 
	 partition function ratios
	 Members:
	 	index
	 	value

*/

struct Pair {
	//pair constructor - nullary
	Pair() {index = value = 0;}

	//define members of pair, index and value
	int index; double value;

	//constructors for non-null
	Pair(int index_, double value_) : index(index_), value(value_) {}
};

//create a cluster structure - stores N points
struct SparseVec {
	//constructor and destructors
	SparseVec(int N) {pairs = new Pair[N];}
	SparseVec() {pairs = NULL;}
	~SparseVec() {delete []pairs;}

	void setNumPoints(int N) {pairs = new Pair[N];}
	Pair* pairs;

	//define bracket operators - allows indexing of points
	Pair& operator[](int index) {return pairs[index];}
	const Pair& operator[](int index) const {return pairs[index];}
};












}