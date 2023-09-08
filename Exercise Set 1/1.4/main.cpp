
#include <chrono>
#include "hoch_rk4.h"


int main (int argc, char * const argv[])
{
    rk4_solution rk4;
    double t0;
    double y0 = 0.0;

    for (int i = 0; i<10; i++)
    {
        rk4_scalar(t0, y0, 0.025, rk4);
        y0 = rk4.y;
        t0 += 0.025;
    }
    

    return 0;
    
}