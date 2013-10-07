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

const int init_pop = 50;
const int no_gens = 100;

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

/* GA parameters */
double mutationRate = 0.015;
int tournamentSize = 5;
bool elitism = true;

class Tour{

    public:
        // Holds our tour of cities
        vector<int> tour;
        int tour_size;

        // Cache
        double fitness;
        double distance;

        Tour(){
            fitness = 0.0;
            distance = 0.0;
            tour_size = n_cities;
            for (int i = 0; i < n_cities; i++) {
                tour.push_back(-1);
            }
        }

        Tour(vector<int> inp_tour){
            fitness = 0.0;
            distance = 0.0;
            tour = inp_tour;
        }

        // Gets a city from the tour
        int getCity(int tourPosition) {
            return tour.at(tourPosition);
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
                for (int cityIndex=0; cityIndex < tour_size; cityIndex++) {
                    // Get city we're travelling from
                    int fromCity = getCity(cityIndex);
                    // City we're travelling to
                    int destinationCity;
                    // Check we're not on our tour's last city, if we are set our
                    // tour's final destination city to our starting city
                    if(cityIndex+1 < tour_size){
            cout<<"ge\n";
                        destinationCity = getCity(cityIndex+1);
                    }
                    else{
                        destinationCity = 0;
                    }
            cout<<cityIndex<<" "<<fromCity<<" "<<destinationCity<<"\n";
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
        string toString() {
            stringstream ret;

            ret<<"|";
            for (int i = 0; i < tour_size; i++) {
                ret<<getCity(i)<<"|";
            }
            return ret.str();
        }
};

class Population {

    public:
        // Holds population of tours
        vector<Tour*> tours;
        int population_size;

        // Construct a population
        Population(int populationSize, bool initialise) {
            population_size = populationSize;
            for (int i = 0; i < population_size; i++) {
                tours.push_back(NULL);
            }
            // If we need to initialise a population of tours do so
            if (initialise) {
                // Loop and create individuals
                for (int i = 0; i < population_size; i++) {
                    Tour* newTour = new Tour();
                    newTour->generateIndividual();
                    saveTour(i, newTour);
                }
            }
        }
        // Saves a tour
        void saveTour(int index, Tour* tour) {
            tours.at(index) = tour;
        }
        // Gets a tour from population
        Tour* getTour(int index) {
            return tours.at(index);
        }

        // Gets the best tour in the population
        Tour* getFittest() {
            Tour* fittest = tours.at(0);
            // Loop through individuals to find fittest
            for (int i = 1; i < population_size; i++) {
                if (fittest->getFitness() <= getTour(i)->getFitness()) {
                    fittest = getTour(i);
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
    int startPos = (int) (r * parent1->tour_size);
    r = ((double)rand()/RAND_MAX);
    int endPos = (int) (r * parent1->tour_size);

    // Loop and add the sub tour from parent1 to our child
    for (int i = 0; i < child->tour_size; i++) {
        // If our start position is less than the end position
        if (startPos < endPos && i > startPos && i < endPos) {
            child->setCity(i, parent1->getCity(i));
        } // If our start position is larger
        else if (startPos > endPos) {
            if (!(i < startPos && i > endPos)) {
                child->setCity(i, parent1->getCity(i));
            }
        }
    }

    // Loop through parent2's city tour
    for (int i = 0; i < parent2->tour_size; i++) {
        // If child doesn't have the city add it
        if (!child->containsCity(parent2->getCity(i))) {
            // Loop to find a spare position in the child's tour
            for (int ii = 0; ii < child->tour_size; ii++) {
                // Spare position found, add city
                if (child->getCity(ii)) {
                }
                else {
                    child->setCity(ii, parent2->getCity(i));
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
    for(int tourPos1=0; tourPos1 < tour->tour_size; tourPos1++){
        // Apply mutation rate
        double r = ((double)rand()/RAND_MAX);
        if(r < mutationRate){
            // Get a second random position in the tour
            r = ((double)rand()/RAND_MAX);
            int tourPos2 = (int) (tour->tour_size * r);

            // Swap them around
            tour->setCity(tourPos2, tourPos1);
            tour->setCity(tourPos1, tourPos2);
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
        tournament.saveTour(i, pop->getTour(randomId));
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
        newPopulation->saveTour(0, pop->getFittest());
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
        newPopulation->saveTour(i, child);
    }

    // Mutate the new population a bit to add some new genetic material
    for (int i = elitismOffset; i < newPopulation->population_size; i++) {
        mutate(newPopulation->getTour(i));
    }

    for(int i=0; i<newPopulation->population_size;i++) {
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

    // Evolve population for 50 generations
    Population* popnew = evolvePopulation(pop);
    delete(pop);
    for (int i = 0; i < no_gens; i++) {
        if(i%2) {
            popnew = evolvePopulation(pop);
            delete(pop);
        }
        else {
            pop = evolvePopulation(popnew);
            delete(popnew);
        }
    }
    cout<<"Final distance: "<<pop->getFittest()->getDistance();
    return 0;
}

