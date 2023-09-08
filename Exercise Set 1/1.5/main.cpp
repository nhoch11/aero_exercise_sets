
#include <chrono>
#include "hoch_rk4_array.h"


int main (int argc, char * const argv[])
{
    rk4_array_solution rk4;
    double t0;
    double* y0 = new double[2];
    // x0
    y0[0] = 0.0;
    y0[1] = 0.0;

    for (int i = 0; i<10; i++)
    {
        rk4_array(t0, y0, 0.025, 2, rk4);
        y0[0] = rk4.y[0];
        y0[1] = rk4.y[1];
        t0 += 0.025;
    }
    

    return 0;
    
}