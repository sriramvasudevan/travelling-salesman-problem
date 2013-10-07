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


/*
 * Global variables which completely specify the input.
 */

bool euclidean ; //true is euclidean.
int n_cities ; //the number of cities.
double **distance_matrix ;
Point *coordinates;

/* GA parameters */
const int init_pop = 100;
const int no_gens = 300;
double mutationRate = 0.015;
int tournamentSize = 15;
bool elitism = true;


class Point {

    public:
        double x;
        double y;
};


class Tour{

    public:
        // Holds our tour of cities
        vector<int> tour;

        // Cache
        double fitness;
        double distance;

        Tour(){
            fitness = 0.0;
            distance = 0.0;
            for (int i = 0; i < n_cities; i++) {
                tour.push_back(-1);
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
            fitness = 0.0;
            distance = 0.0;
        }
        // Gets the tours fitness
        double getFitness() {
            if (fitness == 0.0) {
                fitness = 1/getDistance();
            }
            return fitness;
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

        // Check if the tour contains a city
        bool containsCity(int city){
            if (find(tour.begin(), tour.end(), city) != tour.end()) {
                return true;
            }
            return false;
        }
};


class Population {

    public:
        // Holds population of tours
        vector<Tour*> tours;
        int population_size;

        // Construct a population
        Population(int inp_population_size, bool initialise) {
            population_size = inp_population_size;
            for (int i = 0; i < population_size; i++) {
                tours.push_back(NULL);
            }
            // If we need to initialise a population of tours do so
            if (initialise) {
                // Loop and create individuals
                for (int i = 0; i < population_size; i++) {
                    Tour* newTour = new Tour();
                    newTour->generateIndividual();
                    tours.at(i) = newTour;
                }
            }
        }
        // Gets the best tour in the population
        Tour* getFittest() {
            Tour* fittest = tours.at(0);
            // Loop through individuals to find fittest
            for (int i = 1; i < population_size; i++) {
                if (fittest->getFitness() <= tours.at(i)->getFitness()) {
                    fittest = tours.at(i);
                }
            }
            return fittest;
        }
};


// Applies crossover to a set of parents and creates offspring
Tour* crossover(Tour* parent1, Tour* parent2) {
    // Create new child tour
    Tour* child = new Tour();

    // Get start and end sub tour positions for parent1's tour
    double r = ((double)rand()/RAND_MAX);
    int startPos = (int) (r*n_cities);
    r = ((double)rand()/RAND_MAX);
    int endPos = (int) (r*n_cities);

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
}


// Mutate a tour using swap mutation
void mutate(Tour* tour) {
    // Loop through tour cities
    for(int tourPos1=0; tourPos1 < n_cities; tourPos1++){
        // Apply mutation rate
        double r = ((double)rand()/RAND_MAX);
        if(r < mutationRate){
            // Get a second random position in the tour
            r = ((double)rand()/RAND_MAX);
            int tourPos2 = (int) (n_cities * r);

            // Get the cities at target position in tour
            int city1 = tour->tour.at(tourPos1);
            int city2 = tour->tour.at(tourPos2);

            // Swap them around
            tour->setCity(tourPos2, city1);
            tour->setCity(tourPos1, city2);
        }
    }
}


// Selects candidate tour for crossover
Tour* tournamentSelection(Population* pop) {
    // Create a tournament population
    Population tournament = Population(tournamentSize, false);
    // For each place in the tournament get a random candidate tour and
    // add it
    for (int i = 0; i < tournamentSize; i++) {
        double r = ((double)rand()/RAND_MAX);
        int randomId = (int) (r*pop->population_size);
        tournament.tours.at(i) = pop->tours.at(randomId);
    }
    // Get the fittest tour
    Tour* fittest = tournament.getFittest();
    return fittest;
}


// Evolves a population over one generation
Population* evolvePopulation(Population* pop) {
    Population* newPopulation = new Population(pop->population_size, false);

    // Keep our best individual if elitism is enabled
    int elitismOffset = 0;
    if (elitism) {
        newPopulation->tours.at(0) = pop->getFittest();
        elitismOffset = 1;
    }

    // Crossover population
    // Loop over the new population's size and create individuals from
    // Current population
    for (int i = elitismOffset; i < newPopulation->population_size; i++) {
        // Select parents
        Tour* parent1 = tournamentSelection(pop);
        Tour* parent2 = tournamentSelection(pop);
        // Crossover parents
        Tour* child = crossover(parent1, parent2);
        // Add child to new population
        newPopulation->tours.at(i) =  child;
    }

    // Mutate the new population a bit to add some new genetic material
    for (int i = elitismOffset; i < newPopulation->population_size; i++) {
        mutate(newPopulation->tours.at(i));
    }

    return newPopulation;
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


int main() {
    AcceptInput();
    srand(time(NULL));

    Population* pop = new Population(init_pop, true);
    cout<<"Initial distance: "<<pop->getFittest()->getDistance()<<"\n";

    // Evolve population for no_gens generations
    Population* temp;
    for (int i = 0; i < no_gens; i++) {
        temp = evolvePopulation(pop);
        delete(pop);
        pop = temp;
    }
    Tour* final = pop->getFittest();
    cout<<"Final tour:\n";
    for(int i=0; i<n_cities; i++) {
        cout<<final->tour.at(i)<<" ";
    }
    cout<<"\nFinal distance: "<<final->getDistance()<<"\n";
    return 0;
}
