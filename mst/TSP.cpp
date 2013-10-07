//============================================================================
// Name        : TSP.cpp
// Author      : Viswa
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <time.h>
#include <limits>
#include<cstdlib>

using namespace std;

class Point {

public:
	double x;
	double y;
};

/*
 * Global variables which completely specify the input.
 */
bool euclidean ; //true is euclidean.
int n_cities ; //the number of cities.
double **distance_matrix ;
Point *coordinates;

//The best tour length we've got so far.
double BESTTOURLENGTH = numeric_limits<double>::infinity();

//Given the indices of the tour, it evaluates it.
//If it's better than BESTTOURLENGTH, it updates best tour lengths.
//and outputs the tour to stdout.
void EvalTour(int* tour) {

	double currtourlength = 0 ;
	for ( int i = 1 ; i < n_cities ; i++ ) {
		currtourlength += distance_matrix[tour[i-1]][tour[i]];
	}

	currtourlength += distance_matrix[tour[n_cities-1]][tour[0]] ;
	cout<<"Tour cost is "<<currtourlength<<endl;
	if (currtourlength < BESTTOURLENGTH) {
		BESTTOURLENGTH = currtourlength ;
		//Print tour.
		for ( int i = 0 ; i < n_cities ; i++ ) {
			cout<<(tour[i]+1)<<" ";
		}
		cout<<endl;
	}
}

//Prints the input. For debugging purposes.
void PrintInput() {
	if (!euclidean) cout<<"Not ";
	cout<<"Euclidean"<<endl ;
	cout<<"Number of cities is "<<n_cities<<endl;
	cout<<"Coordinates:"<<endl;
	for ( int i = 0 ; i < n_cities ; i++ ) {
		cout<<"City number "<<i<<": "<<coordinates[i].x<<" "<<coordinates[i].y<<endl;
	}
	cout<<"Distance matrix:"<<endl;
	for ( int i = 0 ; i < n_cities ; i++ ) {
		for ( int j = 0 ; j < n_cities ; j++ ) {
			cout<<distance_matrix[i][j]<<"\t";
		}
		cout<<endl;
	}

}

void AcceptInput() {

	string probtype ;
	cin>>probtype ;
	if (probtype[0]=='e') { //Euclidean.
		euclidean = true ;
	}
	else {
		euclidean = false ;
		cin>>probtype ; //The second word will be euclidean. Scan that off.
	}
	cin>>n_cities ;
	//Put off dynamic allocation.
	distance_matrix = new double* [n_cities] ;
	for ( int i	= 0 ; i < n_cities ; i++ ) {
		distance_matrix[i] = new double[n_cities];
	}
	coordinates = new Point[n_cities] ;

	//Get coordinates.
	for ( int i = 0 ; i < n_cities ; i++ ) {
		cin>>coordinates[i].x>>coordinates[i].y;
	}

	//Fill up distance matrix
	for ( int i = 0 ; i < n_cities ; i++ ) {
		for ( int j = 0 ; j < n_cities ; j++ ) {
			cin>>distance_matrix[i][j] ;
		}
	}
}

// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST
int minKey(double* key, bool* mstSet)
{
   // Initialize min value
   double min = numeric_limits<double>::infinity();

   int min_index;
   for (int i = 0; i < n_cities; i++) {
	  if (mstSet[i] == false && key[i] < min) {
    	 min = key[i];
    	 min_index = i;
     }
   }
   return min_index;
}

// Function to construct and print MST for a graph represented using adjacency
// matrix representation
int* primMST(double** graph)
{
     int* parent = new int[n_cities]; // Array to store constructed MST
     double* key = new double[n_cities];   // Key values used to pick minimum weight edge in cut
     bool* mstSet = new bool[n_cities];  // To represent set of vertices not yet included in MST


     // Initialize all keys as INFINITE
     for (int i = 0; i < n_cities; i++) {
        key[i] = numeric_limits<double>::infinity();
        mstSet[i] = false;
     }

     // Always include a random 1st vertex in MST.

     int randindex = rand()%n_cities ;
     key[randindex] = 0.0;     // Make key 0 so that this vertex is picked as first vertex
     parent[randindex] = -1; // First node is always root of MST
     // The MST will have n_cities vertices
     for (int count = 0; count < n_cities-1; count++)
     {
        // Pick the minimum key vertex from the set of vertices
        // not yet included in MST
    	 int u = minKey(key, mstSet);
        // Add the picked vertex to the MST Set
        mstSet[u] = true;
        // Update key value and parent index of the adjacent vertices of
        // the picked vertex. Consider only those vertices which are not yet
        // included in MST
        for (int v = 0; v < n_cities; v++) {

           // mstSet[v] is false for vertices not yet included in MST
           // Update the key only if graph[u][v] is smaller than key[v]
          if (mstSet[v] == false && graph[u][v] <  key[v]) {
        	  parent[v]  = u;
        	  key[v] = graph[u][v];
          }
        }
     }
     return parent ;
}

void PreorderVisit(int *parent, int *tour, int curr_node) {

	int currtourindex = 0;
	while ( currtourindex<n_cities && tour[currtourindex] != -1) currtourindex++ ;
	tour[currtourindex] = curr_node ;
	//Now to visit its child.
	for ( int i = 0 ; i < n_cities ; i++ ) {
		if (parent[i]==curr_node) PreorderVisit(parent,tour,i);
	}
}

//Using the MST in parent pointer representation,
//outputs the tour. Uses the recursive visit.
void PreorderMSTTour(int *parent, int *tour) {
	int first = 0 ;
	for ( int i = 0 ; i < n_cities ; i++ ) {
		tour[i] = -1 ;

		//Find first node.
		if (parent[i]==-1) first = i ;
	}
	PreorderVisit(parent,tour,first);
}



int main() {
	AcceptInput();
//	PrintInput();
	int *tour = new int[n_cities];
	srand(time(NULL));
	tour[0] = rand()%n_cities ;
	bool flag = true ;

	for ( int i = 1 ; i < n_cities ; ) {
		flag = true ;
		tour[i] = rand()%n_cities ;
		for ( int j = 0 ;  j < i ; j++) {
			if (tour[j]==tour[i]) flag = false ;
		}
		if (flag == true) i++ ;
	}

	EvalTour(tour);
	int *mst_parent_pointer = primMST(distance_matrix);
	PreorderMSTTour(mst_parent_pointer,tour);
	EvalTour(tour);
	return 0;
}
