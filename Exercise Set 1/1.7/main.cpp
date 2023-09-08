
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
 

    // double t00;
    // double* y00 = new double[2];
    // // x0
    // y00[0] = 0.0;
    // y00[1] = 0.0;
    // double* y11 = new double[2];

    // for (int i = 0; i<10; i++)
    // {
    //     rk4_array(t00, y00, 0.025, 2, y11);
        
    //     cout << "t= " << t00+0.025 << "  x= " << setprecision(12) << y11[0] << "  z= " << setprecision(12) << y11[1]<< endl;
    //     y00[0] = y11[0];
    //     y00[1] = y11[1];
    //     t00 += 0.025;

    // }


    return 0;


}
