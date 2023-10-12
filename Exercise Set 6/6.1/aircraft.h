#ifndef aircraft_H
#define aircraft_H


#include "hoch_helper.h"
#include <string>


class aircraft
{
public: // private variables and function can only be accessed inside the class
    aircraft(std::string filename); // class constructor
    ~aircraft(){}; // class destructor (default: releases memory)
    
    void exercise_6_1();

    
private:
    
    // declare a struct from hoch_helper
    Atmosphere m_atm;

    // change all input angles to radians

    // declare json variables
    //simulation
    double m_time_step;
    double m_total_time;

    // aircraft
    double m_wing_area;
    double m_wing_span;
    double m_weight;
    double m_Ixx;
    double m_Iyy;
    double m_Izz;
    double m_Ixy;
    double m_Ixz;
    double m_Iyz;
    double m_hx;
    double m_hy;
    double m_hz;
    double m_cg_shift;
    // put thrust stuff here
    double m_T0;
    double m_T1;
    double m_T2;
    double m_a;

    // initial state
    double m_V;
    double m_altitude;
    double m_elv_angle;
    double m_bank;
    double m_alpha;
    double m_beta;
    double m_p;
    double m_q;
    double m_r;
    double m_heading;
    double m_aileron;
    double m_elevator;
    double m_rudder;
    double m_throttle;

    // aerodynamics
    double m_CL0;
    double m_CL_a;
    double m_CL_qbar;
    double m_CL_de;
    double m_CS_b;
    double m_CS_pbar;
    double m_CS_rbar;
    double m_CS_da;
    double m_CS_dr;
    double m_CDL0;
    double m_CD_L;
    double m_CD_L2;
    double m_CD_S2;
    double m_CD_qbar;
    double m_CD_Lqbar;
    double m_CD_de;
    double m_CD_Lde;
    double m_CD_de2;
    double m_Cl_b;
    double m_Cl_pbar;
    double m_Cl_rbar;
    double m_Cl_Lrbar;
    double m_Cl_da;
    double m_Cl_dr;
    double m_Cm0;
    double m_Cm_a;
    double m_Cm_qbar;
    double m_Cm_de;
    double m_Cn_b;
    double m_Cn_pbar;
    double m_Cn_Lpbar;
    double m_Cn_rbar;
    double m_Cn_da;
    double m_Cn_Lda;
    double m_Cn_dr;

    // additional variables
    double m_cw;
    double m_rho0;

    double m_I[3][3];
    double m_I_inv[3][3];
    
    double m_h[3][3];

    // declare variables for arrow_rk4
    double* m_k1 = new double[12];
    double* m_k2 = new double[12];
    double* m_k3 = new double[12];
    double* m_k4 = new double[12];
    double* m_y_temp = new double[12];


    // functions;
    void aerodynamics_aircraft(double* y, double* ans);

    void aircraft_rk4_func(double t, double* y, double* ans);

    void aircraft_rk4(double t0, double* y0, double dt, int size, double* ans);

    double* init_sim(); // this will initialize the state vector with initial conditions

    

};

#endif