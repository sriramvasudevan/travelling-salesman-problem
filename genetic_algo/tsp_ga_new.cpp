/*
 * tsp_ga.cpp
 * Create a tour and evolve a solution
 */

#include <iostream>
#include <time.h>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <sstream>

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
const int num_nnas = 5 ;

/* GA parameters */
const int init_pop = 100;
const int no_gens = 30;
double mutationRate = 0.15;
int tournamentSize = 15;
int retainedparents = 20;


//The best tour length we've got so far.
double BESTTOURLENGTH = numeric_limits<double>::infinity();

class Tour{

    public:
        // Holds our tour of cities
        vector<int> tour;

        // Cache
        double distance;
        bool updated ; //has the distance been updated?
        Tour(){
            distance = 0.0;
            updated=true;
            for (int i = 0; i < n_cities; i++) {
                tour.push_back(-1);
            }
        }

        Tour(int* i_tour) {
        	distance = 0.0 ;
        	updated=true;
        	for ( int i = 0 ; i < n_cities ; i++ ) {
        		tour.push_back(i_tour[i]);
        	}
        }

        // Creates a random individual
        void generateIndividual() {
            // Loop through all our destination cities and add them to our tour
            for (int cityIndex = 0; cityIndex < n_cities; cityIndex++) {
                setCity(cityIndex, cityIndex);
            }
            // Randomly reorder the tour
            random_shuffle(tour.begin(), tour.end());
        }

        // Sets a city in a certain position within a tour
        void setCity(int tourPosition, int city) {
            tour.at(tourPosition) = city;
            distance = 0.0;
            updated=true;
        }

        // Gets the total distance of the tour
        double getDistance(){
            if (updated) {
                double tourDistance = 0.0;
                // Loop through our tour's cities
                for (int cityIndex=0; cityIndex < n_cities-1; cityIndex++) {
                    tourDistance += distance_matrix[tour[cityIndex]][tour[cityIndex+1]];
                }
                tourDistance += distance_matrix[tour[tour.size()-1]][tour[0]];
                distance = tourDistance;
                updated=false;
            }
            return distance;
        }


        void Eval() {

        	double currtourlength = getDistance();
//        	cout<<currtourlength<<endl;
        	if (currtourlength < BESTTOURLENGTH) {
        		BESTTOURLENGTH = currtourlength ;
        		cout<<currtourlength<<endl ;
        		//Print tour.
        		for ( int i = 0 ; i < n_cities ; i++ ) {
        			cout<<(tour[i]+1)<<" ";
        		}
        		cout<<endl;
        	}

        }

        // Check if the tour contains a city
        bool containsCity(int city){
            if (find(tour.begin(), tour.end(), city) != tour.end()) {
                return true;
            }
            return false;
        }
};

//Function to sort tours.
bool toursortfunction (Tour* a, Tour* b) {
	return (a->getDistance() < b->getDistance()) ;
}


class Population {

    public:
        // Holds population of tours
        vector<Tour*> tours;
        int population_size;

        // Construct a population
        Population(int inp_population_size, vector<Tour*> i_tours ) {

        	Tour* newTour ;
        	for (int i = 0 ; i < i_tours.size() ; i++ ) {
        		newTour = new Tour() ;
        		for (int j = 0 ; j < n_cities ; j++) {
        			newTour->setCity(j,i_tours[i]->tour[j]);
        		}
        		tours.push_back(newTour);
        	}
            population_size = inp_population_size;
            for (int i = i_tours.size(); i < inp_population_size; i++) {
            	newTour = new Tour() ;
            	newTour->generateIndividual();
            	tours.push_back(newTour);
            }
        }

        // Gets the best tour in the population
        Tour* getFittest() {
            Tour* fittest = tours.at(0);
            // Loop through individuals to find fittest
            for (int i = 1; i < tours.size(); i++) {
                if (fittest->getDistance() <= tours.at(i)->getDistance()) {
                    fittest = tours.at(i);
                }
            }
            return fittest;
        }
        vector<Tour*> getFittestFew(int count) {

        	sort (tours.begin(),tours.end(),toursortfunction);
        	vector<Tour*> toReturn (tours.begin(), tours.begin()+count);
        	return toReturn ;

        }
};

// Applies crossover to a set of parents and creates offspring
void crossover(Tour* parent1, Tour* parent2, Tour* child1, Tour* child2) {
    // Get start and end sub tour positions for parent1's tour
    int startPos = rand()%n_cities;
    int endPos = rand()%n_cities;

    if (endPos<startPos) {
        int temp = startPos;
        startPos = endPos;
        endPos = temp;
    }

    int j1=0, j2=0;
    for (int i =0; i<n_cities; i++) {
        if (i>startPos && i<endPos) {
            child1->setCity(j1, parent1->tour.at(i));
            j1++;
            child2->setCity(j2, parent2->tour.at(i));
            j2++;
        }
    }

    int currCityIndex = 0;
    int currCityP1 = 0;
    int currCityP2 = 0;

    for (int i=0; i<n_cities; i++) {
        currCityIndex = (endPos + i)%n_cities;

        currCityP1 = parent1->tour[currCityIndex];
        currCityP2 = parent2->tour[currCityIndex];

        if(!child1->containsCity(currCityP2)) {
            child1->setCity(j1, currCityP2);
            j1++;
        }

        if(!child2->containsCity(currCityP1)) {
            child2->setCity(j2, currCityP1);
            j2++;
        }
    }

    rotate(child1->tour.begin(), child1->tour.end()-startPos-1, child1->tour.end());
    rotate(child2->tour.begin(), child2->tour.end()-startPos-1, child2->tour.end());
}

/*
// Applies crossover to a set of parents and creates offspring
Tour* crossover(Tour* parent1, Tour* parent2) {
    // Create new child tour
    Tour* child = new Tour();

    // Get start and end sub tour positions for parent1's tour
    int startPos = rand()%n_cities;
    int endPos = rand()%n_cities;

    // Loop and add the sub tour from parent1 to our child
    for (int i = 0; i < n_cities; i++) {
        // If our start position is less than the end position
        if (startPos < endPos && i > startPos && i < endPos) {
            child->setCity(i, parent1->tour.at(i));
        }
        // If our start position is larger
        else if (startPos > endPos) {
            if (!(i < startPos && i > endPos)) {
                child->setCity(i, parent1->tour.at(i));
            }
        }
    }

    // Loop through parent2's city tour
    for (int i = 0; i < n_cities; i++) {
    	// If child doesn't have the city add it
        if (!child->containsCity(parent2->tour.at(i))) {
            // Loop to find a spare position in the child's tour
            for (int ii = 0; ii < n_cities; ii++) {
                // Spare position found, add city
                if (child->tour.at(ii)==-1) {
                    child->setCity(ii, parent2->tour.at(i));
                    break;
                }
            }
        }
    }
    return child;
}*/


// Mutate a tour using swap mutation
void mutate(Tour* tour) {
    // Get a random position in the tour
    int tourPos1 = rand()%n_cities;
    // Get a second random position in the tour
    int tourPos2 = rand()%n_cities;
    reverse(tour->tour.begin()+tourPos1, tour->tour.begin()+tourPos2);
    tour->updated = true ;
}


// Selects candidate tour for crossover
Tour* tournamentSelection(Population* pop) {

	//The vector of tours, out of which we'll choose the best one.
	vector<Tour*> tours ;
	int randomId = rand()%(pop->population_size);
	tours.push_back(pop->tours.at(randomId));
	double currbestdistance = tours[0]->getDistance();
	int currbestindex = 0 ;
	for (int i = 1; i < tournamentSize; i++) {
		randomId = rand()%(pop->population_size);
        tours.push_back(pop->tours.at(randomId));
        if (tours[i]->getDistance()<currbestdistance) {
        	currbestdistance=tours[i]->getDistance();
        	currbestindex = i;
        }
    }

	//return the best tour
	return tours[currbestindex];

}


// Evolves a population over one generation
Population* evolvePopulation(Population* pop) {

	vector<Tour*> initialvector = pop->getFittestFew(retainedparents);

	initialvector[0]->Eval();
	Population* newPopulation = new Population(pop->population_size, initialvector);

    // Crossover population
    // Loop over the new population's size and create individuals from
    // Current population
    for (int i = retainedparents; i < newPopulation->population_size; i++) {
        // Select parents
        Tour* parent1 = tournamentSelection(pop);
        Tour* parent2 = tournamentSelection(pop);
        // Create new child tour
        Tour* child1 = new Tour();
        Tour* child2 = new Tour();
        // Crossover parents
        crossover(parent1, parent2, child1, child2);
        // Add child to new population
        newPopulation->tours.at(i) =  child1;
        i++;
        newPopulation->tours.at(i) =  child2;
    }
    for (int i = retainedparents; i < newPopulation->population_size; i++) {
        double r = ((double)rand()/RAND_MAX);
        if(r < mutationRate){
            mutate(newPopulation->tours.at(i));
        }
    }
    vector<Tour*> finaltourvector, newvec ;
    finaltourvector = newPopulation->getFittestFew(pop->population_size/2);
    newvec = pop->getFittestFew(pop->population_size/2);
    finaltourvector.insert(finaltourvector.begin(),newvec.begin(),newvec.end());
    delete(newPopulation);
    Population* toReturn = new Population(pop->population_size,finaltourvector);
    return toReturn;
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

bool Contains (const vector<int> &tour, int tocheck ) {

	if (find(tour.begin(), tour.end(), tocheck) != tour.end()) {
		return true;
	}
	return false;
}


//Get the closest node to tour.size() - 1
//Which is not already in tour.
//Supporter for NNA. Assumes indices etc are proper.
int getClosestNode (const vector<int> &tour) {

	int currcity = tour[tour.size()-1];
	int closestcity = -1 ;
	double leastdistance = numeric_limits<double>::infinity();
	for ( int i = 0 ; i < n_cities ; i++ ) {
		if (!Contains(tour,i) && distance_matrix[currcity][i]<leastdistance) {
			leastdistance = distance_matrix[currcity][i] ;
			closestcity = i ;
		}
	}
	return closestcity ;
}


//Nearest neighbour algorithm.
Tour* NNA () {

	Tour* toReturn = new Tour;
	toReturn->tour.clear();
	toReturn->tour.push_back(rand()%n_cities);
	for ( int i = 1 ; i < n_cities ; i++ ) {
		toReturn->tour.push_back(getClosestNode(toReturn->tour));
	}
	return toReturn ;
}

//Will swapping the edges at these indices help?
bool EdgeExchangeHelps(const vector<int> &it, int e1, int e2 ) {

	//Let's say I have 0 1 2 3 4 5
	//And I want to exchange 0 and 4 (edge node is understood to represent the starting one)
	//New tour would be 0 4 3 2 1 5. I have replaced 0-1 and 4-5 by 0-4 and 1-5
	int s1=e1+1, s2=(e2+1)%(n_cities) ; //successors
	return (distance_matrix[it[e1]][it[e2]] + distance_matrix[it[s1]][it[s2]]+0.1	) < (distance_matrix[it[e1]][it[s1]] + distance_matrix[it[e2]][it[s2]]);

}

//Tries all 4 possible 3 exchanges. Returns the best of the lot..
//Returns 0 if none of them help.
int ThreeEdgeExchangeHelps(const vector<int> &it, int e1, int e2, int e3 ) {

	int toReturn = 0 ;
	int s1=e1+1, s2=(e2+1), s3=(e3+1) ; //successors
	double currbest = distance_matrix[it[e1]][it[s1]] + distance_matrix[it[e2]][it[s2]] + distance_matrix[it[e3]][it[s3]];
	double currcost ;
	//case 1: [e1,e2],[s1,e3],[s2,s3]
	currcost = distance_matrix[it[e1]][it[e2]]+distance_matrix[it[s1]][it[e3]]+distance_matrix[it[s3]][it[s2]];
	if (currcost +0.2< currbest) {
		toReturn = 1 ;
		currbest = currcost ;
	}

	//case2: [e1,s2],[s1,s3],[e2,e3]
	currcost = distance_matrix[it[e1]][it[s2]]+distance_matrix[it[s1]][it[s3]]+distance_matrix[it[e2]][it[e3]];
	if (currcost +0.2 < currbest) {
		toReturn = 2 ;
		currbest = currcost ;
	}

	//case 3: [e1,e3],[s1,s2],[e2,s3]
	currcost = distance_matrix[it[e1]][it[e3]]+distance_matrix[it[s1]][it[s2]]+distance_matrix[it[e2]][it[s3]];
	if (currcost +0.2< currbest) {
		toReturn = 3 ;
		currbest = currcost ;
	}

	//case 4: [e1,s2],[s3,e2],[s1,e3]
	currcost = distance_matrix[it[e1]][it[s2]]+distance_matrix[it[s3]][it[e2]]+distance_matrix[it[s1]][it[e3]];
	if (currcost + 0.2 < currbest) {
		toReturn = 4 ;
		currbest = currcost ;
	}
	return toReturn ;
}

//does the exchange. 4 is from 1-4 (defined by above fn).
void ThreeEdgeExchange (vector<int> &i_tour, int e1, int e2, int e3, int which) {

	vector<int> copy (i_tour); //a copy.
	int s1=e1+1, s2=(e2+1), s3=(e3+1) ; //successors
	switch(which) {

	//case 1: [e1,e2],[s1,e3],[s2,s3]
	case 1:
		i_tour.clear();
		for ( int i = 0 ; i <= e1 ; i++ ) {
			i_tour.push_back(copy[i]);
		}
		for (int i = e2;i>=s1;i--) {
			i_tour.push_back(copy[i]);
		}
		for (int i = e3 ; i >=s2 ; i--) {
			i_tour.push_back(copy[i]);
		}
		for (int i = s3 ; i < n_cities ; i++) {
			i_tour.push_back(copy[i]);
		}
		break;

	//case2: [e1,s2],[s1,s3],[e2,e3]
	case 2:
		i_tour.clear();
		for ( int i = 0 ; i <= e1 ; i++ ) {
			i_tour.push_back(copy[i]);
		}
		for ( int i = s2 ;i<=e3;i++) {
			i_tour.push_back(copy[i]);
		}
		for ( int i = e2 ; i >= s1 ; i--) {
			i_tour.push_back(copy[i]);
		}
		for (int i = s3 ; i < n_cities ; i++) {
			i_tour.push_back(copy[i]);
		}
		break;

	//case 3: [e1,e3],[s1,s2],[e2,s3]
	case 3:
		i_tour.clear();
		for ( int i = 0 ; i <= e1 ; i++ ) {
			i_tour.push_back(copy[i]);
		}
		for (int i = e3 ; i >= s2 ; i--) {
			i_tour.push_back(copy[i]);
		}
		for (int i = s1 ; i <= e2 ; i++) {
			i_tour.push_back(copy[i]);
		}
		for (int i = s3 ; i < n_cities ; i++) {
			i_tour.push_back(copy[i]);
		}
		break;

	//case 4: [e1,s2],[s3,e2],[s1,e3]
	case 4:
		i_tour.clear();
		for ( int i = 0 ; i <= e1 ; i++ ) {
			i_tour.push_back(copy[i]);
		}
		for (int i = s2 ; i <= e3 ; i++) {
			i_tour.push_back(copy[i]);
		}
		for (int i = s1 ; i <= e2 ; i++) {
			i_tour.push_back(copy[i]);
		}
		for (int i = s3 ; i < n_cities ; i++) {
			i_tour.push_back(copy[i]);
		}
		break;
	}
	copy.clear();
}

//Do a two-opt for all the goodtours.
void TwoOpt (vector<Tour*> &goodtours) {

	double prevbest, currbest ;
	vector<double> distvect ;
	bool flag ;
	//for each tour.
	for ( int i = 0 ; i < goodtours.size() ; i++ ) {
		flag = true ;
		while (flag) {
			flag = false ;
			for ( int edge1 = 0 ; edge1 < n_cities -1 ; edge1++ ) {
				for (int edge2 = edge1+2 ; edge2<n_cities ; edge2++ ) {
					if (EdgeExchangeHelps(goodtours[i]->tour,edge1,edge2)) {
						flag = true ;
						//If the tour is 0 1 2 3 4 5 6, and we find that we must exchange when
						//edge1 is 1 and edge2 is 4, our final tour must be 0 1 4 3 2 5 6
						reverse(goodtours[i]->tour.begin()+edge1+1,goodtours[i]->tour.begin()+edge2+1);
						goodtours[i]->updated = true ; //We must set this to true whenever we modify the tour.
						goodtours[i]->Eval();
					}
				}
			}
		}
	}

}
//Do a two-opt for all the goodtours.
void ThreeOpt (vector<Tour*> &goodtours) {

	double prevbest, currbest ;
	vector<double> distvect ;
	bool flag ;
	//for each tour.
	for ( int i = 0 ; i < goodtours.size() ; i++ ) {
		flag = true ;
		//While 3-opt is generating an improvement.
		while (flag) {
			flag = false ;
			for ( int edge1 = 0 ; edge1 < n_cities -1 ; edge1++ ) {
				for (int edge2 = edge1+2 ; edge2<n_cities-1; edge2++ ) {
					for (int edge3 = edge2+2 ; edge3<n_cities-1; edge3++) {
						int which = ThreeEdgeExchangeHelps(goodtours[i]->tour,edge1,edge2,edge3) ;
						if (which>0) {
							flag = true ;
							//If the tour is 0 1 2 3 4 5 6, and we find that we must exchange when
							//edge1 is 1 and edge2 is 4, our final tour must be 0 1 4 3 2 5 6
							ThreeEdgeExchange(goodtours[i]->tour,edge1,edge2,edge3,which);
							goodtours[i]->updated = true ;
							goodtours[i]->Eval();
						}
					}
				}
			}
		}
		cout<<"Optimal 3-opt found for tour no "<<i<<endl;
	}

}


int main() {
    AcceptInput();
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

	Tour(tour).Eval();

	int *mst_parent_pointer = primMST(distance_matrix);
	PreorderMSTTour(mst_parent_pointer,tour);
	Tour *msttour = new Tour(tour);
	msttour->Eval();
	vector<Tour*> initialtours ;
	for (int i = 0 ; i < num_nnas ; i++ ) {
		initialtours.push_back(NNA());
		initialtours[initialtours.size()-1]->Eval();
	}
	initialtours.push_back(msttour);
	TwoOpt(initialtours);
	ThreeOpt(initialtours);
	Population* pop = new Population(init_pop, initialtours);

	// Evolve population for no_gens generations
    Population* temp;
    for (int i = 0; i < no_gens; i++) {
        cout<<"Generation no "<<i<<endl;
    	temp = evolvePopulation(pop);
        delete(pop);
        pop = temp;
    }
    return 0;
}

