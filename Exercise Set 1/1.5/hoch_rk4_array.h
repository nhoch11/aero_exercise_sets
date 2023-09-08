
#ifndef rk4_array_h
#define rk4_array_h
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

//make a struct for the y and t value
struct rk4_array_solution
{
    double* k1;
    double* k2;
    double* k3;
    double* k4;
    double* y;
    int size;
};


void rk4_array(double t0, double* y0, double dt, int size, rk4_array_solution& rk4);

double* f_array(double t, double* y);

#endif