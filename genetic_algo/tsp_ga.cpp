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
const int no_gens = 300;
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

        Tour(){
            distance = 0.0;
            for (int i = 0; i < n_cities; i++) {
                tour.push_back(-1);
            }
        }

        Tour(int* i_tour) {
        	distance = 0.0 ;
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
            // If the tours been altered we need to reset the fitness and distance
            distance = 0.0;
        }

        // Gets the total distance of the tour
        double getDistance(){
            if (distance == 0.0) {
                double tourDistance = 0.0;
                // Loop through our tour's cities
                for (int cityIndex=0; cityIndex < n_cities; cityIndex++) {
                    // Get city we're travelling from
                    int fromCity = tour.at(cityIndex);
                    // City we're travelling to
                    int destinationCity;
                    // Check we're not on our tour's last city, if we are set our
                    // tour's final destination city to our starting city
                    if(cityIndex+1 < n_cities){
                        destinationCity = tour.at(cityIndex+1);
                    }
                    else{
                        destinationCity = 0;
                    }
                    // Get the distance between the two cities
                    tourDistance += distance_matrix[fromCity][destinationCity];
                }
                distance = tourDistance;
            }
            return distance;
        }

        void Eval() {

        	double currtourlength = getDistance();
        	cout<<currtourlength<<endl;
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

// Applies crossover to a set of parents and creates offspring
Tour* crossover(Tour* parent1, Tour* parent2) {
    // Create new child tour
    Tour* child1 = new Tour();
    Tour* child2 = new Tour();

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

    return child1;
}


// Mutate a tour using swap mutation
void mutate(Tour* tour) {
    // Get a random position in the tour
    int tourPos1 = rand()%n_cities;
    // Get a second random position in the tour
    int tourPos2 = rand()%n_cities;
    reverse(tour->tour.begin()+tourPos1, tour->tour.begin()+tourPos2);
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
    newPopulation->tours.at(0) = pop->getFittest();

    // Crossover population
    // Loop over the new population's size and create individuals from
    // Current population
    for (int i = retainedparents; i < newPopulation->population_size; i++) {
        // Select parents
        Tour* parent1 = tournamentSelection(pop);
        Tour* parent2 = tournamentSelection(pop);
        // Crossover parents
        Tour* child = crossover(parent1, parent2);
        // Add child to new population
        newPopulation->tours.at(i) =  child;
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

	for (int i = 0 ; i < tour.size();i++) {
		if (tour[i]==tocheck) return true ;
	}
	return false ;
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

//Do a two-opt for all the goodtours.
void TwoOpt (vector<Tour*> &goodtours) {

	double prevbest, currbest ;
	//for each tour.
	for ( int i = 0 ; i < goodtours.size() ; i++ ) {
		do {
			prevbest = goodtours[i]->getDistance();
			for ( int city = 0 ; city < n_cities ; city ++ ) {

			}

		} while(currbest < prevbest );

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
	vector<Tour*> initialtoursforga ;
	for (int i = 0 ; i < num_nnas ; i++ ) {
		initialtoursforga.push_back(NNA());
		initialtoursforga[initialtoursforga.size()-1]->Eval();
	}
	initialtoursforga.push_back(msttour);
	Population* pop = new Population(init_pop, initialtoursforga);

	// Evolve population for no_gens generations
    Population* temp;
    for (int i = 0; i < no_gens; i++) {
        temp = evolvePopulation(pop);
        delete(pop);
        pop = temp;
    }
    vector<Tour*> toursforopt = pop->getFittestFew(10);
    TwoOpt(toursforopt);
    return 0;
}

