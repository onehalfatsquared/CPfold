#include <cstdlib>
#include <stdio.h>
#include "bDynamics.h"



int main(int argc, char* argv[]) {

	//handle input
	if (argc != 3) {
		fprintf(stderr, "Usage: %s <N> <state>", argv[0]);
	}
	int N = atoi(argv[1]);
	int state = atoi(argv[2]);

	//get the database here? pass to estimate mfpt?














}