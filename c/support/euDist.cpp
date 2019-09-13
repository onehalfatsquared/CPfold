#include <math.h>
#include "adjacency.h"
#include "../defines.h"
namespace bd {

void c2p(double* cluster, double* particles, int N) {
	//converts array of form (x1,y1,x2,y2,...) to (x1,x,2,x3,x4,...)

#if (DIMENSION == 2) 
	for (int i = 0; i < 2*N; i++) {
		if (i % 2 == 0) {
			particles[i/2] = cluster[i];
		}
		else {
			particles[N+i/2] = cluster[i];
		}
	}
#elif (DIMENSION == 3)
	for (int i = 0; i < 3*N; i++) {
		if (i % 3 == 0) {
			particles[i/3] = cluster[i];
		}
		else if (i % 3 == 1) {
			particles[N+i/3] = cluster[i];
		}
		else {
			particles[2*N+i/3] = cluster[i];
		}
	}
#endif

}


double euDist(double* particles, int i, int j, int N, double* Z) {
	//computes euclidean distance between particles i and j
  Z[0] = particles[i]-particles[j];
  Z[1] = particles[N+i]-particles[N+j];
#if (DIMENSION == 3)
	Z[2] = particles[2*N+i]-particles[2*N+j];
#endif


  double R = Z[0]*Z[0]+Z[1]*Z[1];
#if (DIMENSION == 3)
  R += Z[2]*Z[2];
#endif

  R = sqrt(R);
  Z[0] /= R; Z[1] /= R;
#if (DIMENSION == 3)
  Z[2] /= R;
#endif

  return R;
}


}
