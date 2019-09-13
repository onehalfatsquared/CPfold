#include <math.h>
#include <algorithm>
#include "bDynamics.h"
namespace bd {

double ljP(double r, double rho, double E){
  //evaluate the lennard-jones potential derivative
  double Y = pow(1/r, rho);
  return -4 * E * rho / r * Y * (2*Y-1);
}

double ljEval(double* particles, int rho, double* E, int N, int* P) {
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
                    S += ljP(r, rho, Eeff); 
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

void ljGrad(double* particles, int rho, double* E, int N, int* P, double* g) {
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
            S[k] = S[k] + ljP(r, rho, Eeff) * Z[k];
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