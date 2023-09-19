#ifndef ARROW_H
#define ARROW_H


#include "hoch_helper.h"
#include <string>


class arrow
{
public: // private variables and function can only be accessed inside the class
    arrow(std::string filename); // class constructor
    ~arrow(){}; // class destructor (default: releases memory)

    void exercise_2_2();
    void exercise_2_4();
    
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
    double m_Iyy;
    // aerodynamic coefficients
    double m_CLa;
    double m_CD0;
    double m_CD2;
    double m_Cma;
    double m_Cmq;

    // // declare calculated variables
    // double m_Fxb;
    // double m_Fzb;
    // double m_Myb;
    // double m_CL;
    // double m_CD;
    // double m_Cm;
    
    // double m_alpha;
    // double m_V;
    // double m_g;
    // double m_rho;
    
    // declare variables for arrow_rk4
    double* m_k1 = new double[12];
    double* m_k2 = new double[12];
    double* m_k3 = new double[12];
    double* m_k4 = new double[12];
    double* m_y_temp = new double[12];


    // declare arrays for holding info
    // double* m_y;
    // double* m_dy;

    // functions;
    void arrow_rk4_2_2(double t0, double* y0, double dt, int size, double* ans);
    void arrow_rk4_func_2_2(double t, double*y, double* ans);
    void aerodynamics_2_2(double* y, double* ans);
    
    void arrow_rk4_2_4(double t0, double* y0, double dt, int size, double* ans);
    void arrow_rk4_func_2_4(double t, double*y, double* ans);
    void aerodynamics_2_4(double* y, double* ans);
    

};

#endif