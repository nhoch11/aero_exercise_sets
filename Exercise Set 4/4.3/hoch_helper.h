
#ifndef atmospheric_properties_h
#define atmospheric_properties_h
#include "cmath"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "json.hpp"

using namespace std; // use standard or 'std' namespace
using json = nlohmann::json;

// Define a constant for pi
#define pi 3.141592653589793238462643383279


struct Atmosphere // Atmosphere container has these things in it
{
    double geopotential_altitude;
    double temperature;
    double pressure;
    double density;
    double speed_of_sound;
    double dynamic_viscosity;
};


void get_atmospheric_properties_si(double altitude, Atmosphere& atm);

void get_atmospheric_properties_english(double altitude, Atmosphere& atm);

void rk4_array(double t0, double* y0, double dt, int size, double* ans);

void array_copy(double* A, double* B, int size);

void quat_mult(double* A, double* B, double* ans);

void quat_norm(double* quat);

void euler_to_quat(double* eul, double* quat);

void quat_to_euler(double* quat, double* eul);

void body_to_fixed(double* vec, double* quat, double* ans);

void fixed_to_body(double* vec, double* quat, double* ans);

double* f_array(double t, double* y);

double gravity_si(double H);

double gravity_english(double H);

#endif