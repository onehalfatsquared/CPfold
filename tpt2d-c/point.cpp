#include "point.h"
namespace bd{

//functions to go to/from arrays and clusters in 2 and 3 dimensions
	
void Cluster::makeArray2d(double* array, int N) {
	for (int j = 0; j < N; j++) {
		array[2*j] = points[j].x;
		array[2*j+1] = points[j].y;
	}
}

void Cluster::makeArray3d(double* array, int N) {
	for (int j = 0; j < N; j++) {
		array[3*j] = points[j].x;
		array[3*j+1] = points[j].y;
		array[3*j+2] = points[j].z;
	}
}

void Cluster::makeCluster2d(double* array, int N) {
	for (int j = 0; j < N; j+=2) {
		double x = array[j]; double y = array[j+1]; 
		points[j] = Point(x,y);
	}
}

void Cluster::makeCluster3d(double* array, int N) {
	for (int j = 0; j < N; j+=3) {
		double x = array[j]; double y = array[j+1]; double z = array[j+2];
		points[j] = Point(x,y,z);
	}
}

}