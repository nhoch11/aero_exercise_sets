#include <iostream>
#include <cmath>
#include <cstdio>

double* get_atmospheric_properties_English(double altitude)
{
    double* ans; 
    ans = new double[5];
    // convert input of feet to meters
    altitude = altitude*0.3048;
    
    // declare constants and initialize variables
    const double RE = 6356766.0; // meters to feet
    const double G0 = 9.806645; // m^2 /s
    double Z = (RE*altitude)/(RE + altitude);
    double p; // pressure Pa
    double T; // Temperature K
    double rho = 0.0; // kg/m^3
    double a = 0.0; // m/s
    // standard sea level conditions
    
    const double R = 287.0528; // Nm/kgK
    
    // store pre calculated temperatures and pressures at each layer boundary
    double T_array[7] = {288.150, 216.650, 216.650, 228.650, 270.650, 270.650, 252.650}; //K
    double p_array[7] = {101325.0, 22632.0491189944,5474.88167406511, 868.016875642432, 110.905974487888,59.0007483456162, 18.2100049756603 }; // Pa
    double T_prime[7] = {-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.002, -0.004}; // k/m
    double Z_const[7] = {0.0, 11000.0, 20000.0, 32000.0, 47000.0, 52000.0, 61000.0}; // m
   

    if (Z <= 11000.0){
        T = T_array[0] + T_prime[0]*(Z-Z_const[0]); // K
        p = p_array[0]*pow(((T_array[0] +T_prime[0]*(Z-Z_const[0]))/T_array[0]),(-G0/(R*T_prime[0])));
    } else if (Z > 11000.0 && Z <= 20000.0){
       // add the portion of the second altitude range
        T = T_array[1] + T_prime[1]*(Z-Z_const[1]); // K
        // use the Ti prime = 0 equation
        p = p_array[1]*exp(-(G0*(Z-Z_const[1]))/(R*T_array[1])); // Pa
        
    } else if (32000.0 >= Z && Z > 20000.0){
        // add the portion of the third altitude range
        T = T_array[2] + T_prime[2]*(Z-Z_const[2]); // K
        // use the Ti prime /= 0 equation
        p = p_array[2]*pow(((T_array[2] +T_prime[2]*(Z-Z_const[2]))/T_array[2]),(-G0/(R*T_prime[2])));
        
    } else if (47000.0 >= Z && Z > 32000.0){
        // add the portion of the third altitude range
        T = T_array[3] + T_prime[3]*(Z-Z_const[3]); // K
        // use the Ti prime /= 0 equation
        p = p_array[3]*pow(((T_array[3] + T_prime[3]*(Z-Z_const[3]))/T_array[3]),(-G0/(R*T_prime[3])));
        
    } else if (52000.0 >= Z && Z > 47000.0){
        // add the portion of the third altitude range
        T = T_array[4] + T_prime[4]*(Z-Z_const[4]); // K
        // use the Ti prime /= 0 equation
        p = p_array[4]*exp((-G0*(Z-Z_const[4]))/(R*T_array[4]));
        
    } else if (61000.0 >= Z && Z > 52000.0){
        // add the portion of the third altitude range
        T = T_array[5] + T_prime[5]*(Z-Z_const[5]); // K
        // use the Ti prime /= 0 equation
        p = p_array[5]*pow(((T_array[5]  + T_prime[5]*(Z-Z_const[5]))/T_array[5] ),(-G0/(R*T_prime[5])));
        
    } else if (79000.0 >= Z && Z > 61000.0){
        // add the portion of the third altitude range
        T = T_array[6] + T_prime[6]*(Z-Z_const[6]); // K
        // use the Ti prime /= 0 equation
        p = p_array[6]*pow(((T_array[6] + T_prime[6]*(Z-Z_const[6]))/ T_array[6]),(-G0/(R*T_prime[6])));
        
    }
    rho = p/(R*T);
    a = sqrt(1.4*R*T);

    ans[0] = Z/0.3048; // meters to ft
    ans[1] = T*1.8;  // Kelvin to Rankine
    ans[2] = p*0.020885434304801722; // Pa to lbf/ft^2
    ans[3] = rho*0.00194032032363104; // kg/m^3 to slugs/ ft^3
    ans[4] = a/0.3048; // m/s to ft/s
    return ans;
}


int atmosphere_print_English()
{
    double* ans;
    FILE* en_file = fopen("stdatmos_English.txt", "w");

    fprintf(en_file, " Geometric_Alt[ft]    Geopotential_Alt[ft]   Temperature[R]     Pressure[lbf/ft^2]    Density[slugs/ft^3]    Speed of Sound[ft/s]\n");

    for (double H = 0.0; H < 180000.00; H += 5000.0) {
        ans = get_atmospheric_properties_English(H);
        fprintf(en_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", H,ans[0], ans[1], ans[2], ans[3], ans[4]);
    }
    fclose(en_file);
    return 0;
}


int main() {
    atmosphere_print_English();
    std::cout << "check" << std::endl;
    return 0;
}