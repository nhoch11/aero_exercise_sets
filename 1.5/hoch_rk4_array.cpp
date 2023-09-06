#include "hoch_rk4_array.h"



void rk4_array(double t0, double* y0, double dt, int size, rk4_array_solution& rk4)
{
    rk4.size = size;
    rk4.k1 = new double[size];
    rk4.k2 = new double[size];
    rk4.k3 = new double[size];
    rk4.k4 = new double[size];

    // calculate ks
    double* y_temp = new double[rk4.size];
    

    for (int i = 0; i < rk4.size; i++)
    {
        double* k1 = f_array(t0, y0);
        
        rk4.k1[i] = dt * k1[i];
    }

    for (int i = 0; i < rk4.size; i++)
    {
        y_temp[i] = y0[i] + 0.5*rk4.k1[i];
    }
        double* k2 = f_array(t0 + 0.5*dt, y_temp);
    
    for (int i = 0; i < rk4.size; i++)
    {
        rk4.k2[i] = dt * k2[i];
    }
        
    for (int i = 0; i < rk4.size; i++)
    {
        y_temp[i] = y0[i] + 0.5*rk4.k2[i];
        
    }
        double* k3 = f_array(t0 + 0.5*dt, y_temp);
        
    for (int i = 0; i < rk4.size; i++)
    {  
        rk4.k3[i] = dt * k3[i];
    }
    
    for (int i = 0; i < rk4.size; i++)
    {
        y_temp[i] = y0[i] + rk4.k3[i];
    }
        double* k4 = f_array(t0 + dt, y_temp);
        
    for (int i = 0; i < rk4.size; i++)
    {
        rk4.k4[i] = dt * k4[i];
    }
    for (int i = 0; i < rk4.size; i++)
    {
        rk4.y[i] = y0[i] + (rk4.k1[i] + 2*rk4.k2[i] + 2*rk4.k3[i] + rk4.k4[i])/6.0  ;
    }
    
    
    cout << "t= " << t0+dt << "  x= " << setprecision(12) << rk4.y[0] << "  z= " << setprecision(12) << rk4.y[1]<< endl;
}

double* f_array(double t, double* y)
{
    double* derivatives = new double[2];

    // x = y[0]     z = y[1]
    // dxdt
    derivatives[0] =  1.0 + y[1]*y[1]*sin(y[0]);
    // dzdt
    derivatives[1] =  1.0 + y[0]*cos(y[1]);
    return derivatives;
}

