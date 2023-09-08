
#include <chrono>
#include "hoch_helper.h"


int main (int argc, char * const argv[])
{
    Atmosphere atm;
    double H;

    auto start_time = chrono::high_resolution_clock::now();
    for (H = 0.0; H < 75000.0; H += 0.005)
    {
        get_atmospheric_properties_si(H,atm);
    }
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end_time - start_time;
    
    printf("SI units total time (sec): %f\n",duration.count());
    
    start_time = chrono::high_resolution_clock::now();
    for (H = 0.0; H < 200000.0; H += 0.005)
    {
        get_atmospheric_properties_english(H,atm);
    }
    end_time = chrono::high_resolution_clock::now();
    duration = end_time - start_time;
    
    printf("English units total time (sec): %f\n",duration.count());
    
    return 0;
    
}