#ifndef CANNONBALL_H
#define CANNONBALL_H


#include "hoch_helper.h"
#include <string>


class cannonball
{
public: // private variables and function can only be accessed inside the class
    cannonball(std::string filename); // class constructor
    ~cannonball(){}; // class destructor (default: releases memory)
    
    void exercise_4_3();
    
private:
    
    // declare a struct from hoch_helper
    Atmosphere m_atm;

    // declare json variables
    double m_time_step;
    double m_ref_area;
    double m_ref_length;
    double m_init_V;
    double m_init_altitude;
    double m_init_theta;
    double m_weight;
    double m_Ixx;
    double m_Iyy;
    // aerodynamic coefficients
    double m_CLa;
    double m_CD0;
    double m_CD2;
    double m_Cma;
    double m_Cmq;
    double m_Clp;
    double m_Cl0;

    
    // declare variables for arrow_rk4
    double* m_k1 = new double[12];
    double* m_k2 = new double[12];
    double* m_k3 = new double[12];
    double* m_k4 = new double[12];
    double* m_y_temp = new double[12];


    // functions;

    void aerodynamics_cannonball(double* y, double* ans);
    double get_sphere_CD(double reynolds);

    void cannonball_rk4_func(double t, double* y, double* ans);
    void cannonball_rk4(double t0, double* y0, double dt, int size, double* ans);

    

};

#endif