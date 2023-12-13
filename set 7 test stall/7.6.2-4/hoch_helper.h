
#ifndef atmospheric_properties_h
#define atmospheric_properties_h
#include "cmath"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "json.hpp"
#include "LUDv.h"
#include <stdio.h>

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

double* cross(double x1, double y1, double z1, double x2, double y2, double z2);

bool matrix_invert_3x3(double a[3][3], double b[3][3]);

void matrix_vector_mult_3(double rm[3][3], double v[3], double* ans);

void array_print(double* arr, int size);

void array_print_3x3(double arr[3][3]);

void vector_cross_3(double a[3], double b[3], double ans[3]);

void vector_normalize_3(double vec[3]);

void matrix_AxB_solve(double** A, double* B, int size, double* x);

#endif