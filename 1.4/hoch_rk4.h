
#ifndef rk4_h
#define rk4_h
#include "cmath"
#include <iostream>
#include <iomanip>

using namespace std;

//make a struct for the y and t value
struct rk4_solution
{
    double k1;
    double k2;
    double k3;
    double k4;
    double y;
};


void rk4_scalar(double t0, double y0, double dt, rk4_solution& rk4);

double f_scalar(double t, double y);

#endif