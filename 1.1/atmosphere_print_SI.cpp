#include <cstdio>
#include <iostream>
#include "get_atmospheric_properties_SI.cpp"

int atmosphere_print_SI()
{
    double* ans;
    FILE* si_file = fopen("stdatmos_si.txt", "w");

    fprintf(si_file, " Geometric_Alt[m]      Geopotential_Alt[m]     Temperature[K]      Pressure[N/m^2]      Density[kg/m^3]      Speed of Sound[m/s]\n");

    for (double H = 0.0; H < 72000.00; H += 2000.0) {
        ans = get_atmospheric_properties_SI(H);
        fprintf(si_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", H,ans[0], ans[1], ans[2], ans[3], ans[4]);
    }
    fclose(si_file);
    std::cout << "check" << std::endl;
    return 0;
}