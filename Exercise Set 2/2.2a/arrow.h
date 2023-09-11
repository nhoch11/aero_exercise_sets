#ifndef ARROW_H
#define ARROW_H


#include "hoch_helper.h"

class arrow
{
public: // private variables and function can only be accessed inside the class
    arrow(string filename); // class constructor
    ~arrow(){}; // class destructor (default: releases memory)


    
private:

    // declare json variables
    double m_time_step;
    double m_ref_area;
    double m_ref_length;
    double m_V;
    double m_z;
    double m_theta;
    double m_weight;
    double m_Iyy;
    double m_CL_a;
    double m_CD0;
    double m_CD2;
    double m_Cm_a;
    double m_Cm_q;

    // declare calculated variables
    double m_F_xb;
    double m_F_zb;
    double m_M_yb;
    double m_CL;
    double m_CD;
    double m_Cm;
    double m_V;
    double m_alpha;
    double m_mass;
    double m_g;

    // declare state variables
    double m_u;
    double m_w;
    double m_q;
    double m_x_f;
    double m_z_f;
    double m_theta;

    // declare arrays for holding info
    double m_y[6];
    double m_dy[6];
    //double* FM[3];

    // functions
    void aerodynamics_2_2(double* y, double* FM);
    
    double* arrow_EOM(double t, double* y);

};

#endif