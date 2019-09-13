#include <math.h>
#include <algorithm>
#include "bDynamics.h"
namespace bd {

double morseP(double r, double rho, double E){
  //evaluate the morse potential derivative
  double Y = exp(-(rho) * ((r) - 1));
  return -2 * (rho) * (E) * ( Y*Y -Y);
}

double morseEval(double* particles, int rho, double* E, int N, int* P) {
  //compute total energy of system with pairwise morse potential //not necc?
	double S = 0; int i, j; int rep = 0;
    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            if(j!=i){
                double* Z = new double[2];
                double r = euDist(particles, i, j, N, Z);
                int minimum = std::min(i,j);
                int maximum = std::max(i,j);
                if(P[N*maximum+minimum]==1){
                    double Eeff = E[N*maximum+minimum];
                    S += morseP(r, rho, Eeff); 
                }
                else if(r<1){
                    S -= rep/r; 
                }
              delete []Z;
            }
        }
    }
    return S;
}

void morseGrad(double* particles, int rho, double* E, int N, int* P, double* g) {
  //compute gradient of energy of system with pairwise morse potential
	int i, j; int rep = 0;
	double* S = new double[2];
	for(i=0;i<N;i++){
    S[0] = S[1] = 0;
    for(j=0;j<N;j++){
      if(j!=i){
        double* Z = new double[2];
        double r = euDist(particles, i, j, N, Z);
        int minimum = std::min(i,j);
        int maximum = std::max(i,j);
        if(P[N*maximum+minimum]==1){
          double Eeff = E[N*maximum+minimum];
          for(int k=0; k<2;k++)
            S[k] = S[k] + morseP(r, rho, Eeff) * Z[k];
        }
        else if(r<1){
          for(int k=0; k<2;k++)
            S[k] = S[k] - rep * Z[k] / (r*r);
        }
      	delete []Z;
    	}
		}
    for(int k=0; k<2;k++)
      g[2*i+k] = S[k];
  }
  delete []S;
}

}