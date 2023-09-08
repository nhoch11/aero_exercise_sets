
#ifndef atmospheric_properties_h
#define atmospheric_properties_h
#include "cmath"
#include <iostream>
#include <iomanip>

using namespace std; // use standard or 'std' namespace

struct Atmosphere // Atmosphere container has these things in it
{
    double geopotential_altitude;
    double temperature;
    double pressure;
    double density;
    double speed_of_sound;
};


void get_atmospheric_properties_si(double altitude, Atmosphere& atm);

void get_atmospheric_properties_english(double altitude, Atmosphere& atm);

void rk4_array(double t0, double* y0, double dt, int size, double* ans);

double* f_array(double t, double* y);

#endif