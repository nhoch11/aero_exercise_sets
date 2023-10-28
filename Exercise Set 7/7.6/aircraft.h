#ifndef aircraft_H
#define aircraft_H


#include "hoch_helper.h"
#include <string>


class aircraft
{
public: // private variables and function can only be accessed inside the class
    aircraft(std::string filename); // class constructor
    ~aircraft(){}; // class destructor (default: releases memory)
    
    void init_sim(); // this will initialize the state vector with initial conditions
    void run_sim();

    
private:
    
    // declare a struct from hoch_helper
    Atmosphere m_atm;

    double* m_initial_state;
    //double* m_trim_state;
    int m_size;
    FILE* m_check_file;

    // declare json variables
    json m_input;

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
    double m_bank_angle;
    double m_alpha;
    double m_beta;
    double m_p;
    double m_q;
    double m_r;
    double m_heading;
    double m_da;
    double m_de;
    double m_dr;
    double m_throttle;

    // initial trim
    string m_trim_type;
    double m_climb_angle;
    double m_sideslip_angle;

    double m_finite_diff_step;
    double m_relaxation;
    double m_tol;
    int m_max_iter;
    bool m_verbose;

    // aerodynamics
    double m_CL0;
    double m_CL_a;
    double m_CL_qbar;
    double m_CL_de;
    double m_CS_b;
    double m_CS_pbar;
    double m_CS_Lpbar;
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

    // create controls vector
    double* m_controls;

    double m_thrust_mag;

    double m_thrust_loc[3];


    double m_thrust_dir[3];

    double m_CG_shift[3];
    


    // functions;

    void init_from_state();

    void init_from_trim();

    void calc_R(double G[6], double* y, double phi, double theta, double R[6]);

    void aerodynamics_aircraft(double* y, double* ans);

    void aircraft_rk4_func(double t, double* y, double* ans);

    void aircraft_rk4(double t0, double* y0, double dt, int size, double* ans);

    
    

};

#endif