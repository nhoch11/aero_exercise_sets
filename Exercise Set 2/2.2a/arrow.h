#ifndef ARROW_H
#define ARROW_H

#include "hoch_helper.h"

class arrow
{
public: // private variables and function can only be accessed inside the class
    arrow(string filename); // class constructor


    


    // functions
    double* m_arrow_eqs(double t, double* y);

    ~arrow(){}; // class destructor (default: releases memory)

private:

    // variables

    double m_time_step;
    double m_u;
    double m_w;
    double m_q;
    double m_x_f;
    double m_z_f;
    double m_theta;

};

#endif