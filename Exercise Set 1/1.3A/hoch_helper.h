
#ifndef atmospheric_properties_h
#define atmospheric_properties_h
#include "cmath"
#include <iostream>

using namespace std;
// make a containter called atmosphere (at a specific memory locaiton)
struct Atmosphere
// Atmosphere container has these things in it
{
    double geopotential_altitude;
    double temperature;
    double pressure;
    double density;
    double speed_of_sound;
};





void get_atmospheric_properties_si(double altitude, Atmosphere& atm);

void get_atmospheric_properties_english(double altitude, Atmosphere& atm);


#endif