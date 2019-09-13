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














}