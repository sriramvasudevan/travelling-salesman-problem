/*
 * City.cpp
 * Models a city
 */
#include <math.h>

class City {
    int x;
    int y;

    public:
    // Constructs a randomly placed city
    City(){
        this.x = (int)(Math.random()*200);
        this.y = (int)(Math.random()*200);
    }
    // Constructs a city at chosen x, y location
    City(int x, int y){
        this.x = x;
        this.y = y;
    }
    // Gets city's x coordinate
    int getX(){
        return x;
    }
    // Gets city's y coordinate
    int getY(){
        return y;
    }
    // Gets the distance to given city
    double distanceTo(City city){
        int xDistance = abs(getX() - city.getX());
        int yDistance = abs(getY() - city.getY());
        double distance = sqrt( (xDistance*xDistance) + (yDistance*yDistance) );
        return distance;
    }
    string toString(){
        return getX()+", "+getY();
    }
}


