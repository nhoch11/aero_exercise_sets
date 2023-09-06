#include <iostream>
#include <chrono>
#include "hoch_helper.h"


int main (int argc, char * const argv[])
{
    double* ans;
    double H;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    for (H = 0.0; H < 75000.0; H += 0.005)
    {
        ans = get_atmospheric_properties_si(H);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    
    printf("SI units total time (sec): %f\n",duration.count());
    
    start_time = std::chrono::high_resolution_clock::now();
    for (H = 0.0; H < 200000.0; H += 0.005)
    {
        ans = get_atmospheric_properties_english(H);
    }
    end_time = std::chrono::high_resolution_clock::now();
    duration = end_time - start_time;
    
    printf("English units total time (sec): %f\n",duration.count());
    
    return 0;
    
}