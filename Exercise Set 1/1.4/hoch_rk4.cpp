#include "hoch_rk4.h"



void rk4_scalar(double t0, double y0, double dt, rk4_solution& rk4)
{

    // calculate ks
    rk4.k1 = dt * f_scalar(t0, y0);
    rk4.k2 = dt * f_scalar(t0 + 0.5*dt, y0 + 0.5*rk4.k1);
    rk4.k3 = dt * f_scalar(t0 + 0.5*dt, y0 + 0.5*rk4.k2);
    rk4.k4 = dt * f_scalar(t0 + dt, y0 + rk4.k3);
    
    rk4.y = y0 + (rk4.k1 + 2*rk4.k2 + 2*rk4.k3 + rk4.k4)/6.0  ;  
    
    cout << "t= " << setprecision(15) << t0+dt << "  y= " << setprecision(15) << rk4.y << endl;
}

double f_scalar(double t, double y)
{
    return 1 + tan(y);
}

