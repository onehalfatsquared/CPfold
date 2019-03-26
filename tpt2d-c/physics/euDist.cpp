#include <math.h>
#include "bDynamics.h"
namespace bd {

void c2p(double* cluster, double* particles, int N) {
	//converts array of form (x1,y1,x2,y2,...) to (x1,x,2,x3,x4,...)
	for (int i = 0; i < 2*N; i++) {
		if (i % 2 == 0) {
			particles[i/2] = cluster[i];
		}
		else{
			particles[N+i/2] = cluster[i];
		}
	}
}

double euDist(double* particles, int i, int j, int N, double* Z){
	//computes euclidean distance between particles i and j
    Z[0] = particles[i]-particles[j];
    Z[1] = particles[N+i]-particles[N+j];

    double R = sqrt(Z[0]*Z[0]+Z[1]*Z[1]);
    Z[0] /= R; Z[1] /= R;
    return R;
}

}