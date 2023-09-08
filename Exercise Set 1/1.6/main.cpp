
#include <chrono>
#include "hoch_helper.h"

int main (int argc, char * const argv[])
{
    double t=0;
    double dt = 0.0001;
    double* y = new double[2];
    double* y1 = new double[2];
    y[0] = 0.0;
    y[1] = 0.0;

    auto start_time = chrono::high_resolution_clock::now();
    for (int i = 0; i < 10000000; i++)
    {
        rk4_array(t,y,dt,2,y1);
        y[0] = y1[0];
        y[1] = y1[1];
        t += dt;
    }
    auto end_time = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    printf("RK4 array total time (sec): %f\n",duration.count());
 
    return 0;
}
