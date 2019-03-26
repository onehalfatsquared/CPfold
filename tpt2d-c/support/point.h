#pragma once
#include <cstdlib>

namespace bd {

//create a point structure - stores x, y, z values of a point
struct Point {
	//point constructor - nullary
	Point() {x = y = z = 0;}

	//define members of point, coordinates
	double x, y, z;

	//constructors for non-null
	Point(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
	Point(double x_, double y_) : x(x_), y(y_), z(0) {}
};

//create a cluster structure - stores N points
struct Cluster {
	//constructor and destructors
	Cluster(int N) {points = new Point[N];}
	Cluster() {points = NULL;}
	~Cluster() {delete []points;}

	void setNumPoints(int N) {points = new Point[N];}
	Point* points;

	//define bracket operators - allows indexing of points
	Point& operator[](int index) {return points[index];}
	const Point& operator[](int index) const {return points[index];}

	//functions to go from/to arrays and clusters
	void makeArray2d(double* array, int N); 
	void makeCluster2d(double* array, int N);
	void makeArray3d(double* array, int N); 
	void makeCluster3d(double* array, int N);
};

}