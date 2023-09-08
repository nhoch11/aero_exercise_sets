#include <iostream>
#include <cmath>

double* get_atmospheric_properties_SI(double altitude)
{
    double* ans; 
    ans = new double[5];
    // declare constants and initialize variables
    const double RE = 6356766.0; // meters
    const double G0 = 9.806645; // m^2 /s
    double Z = (RE*altitude)/(RE + altitude);
    double T = 0.0; //K
    double p = 0.0; // Pa
    double rho = 0.0; // kg/m^3
    double a = 0.0; // m/s
    // standard sea level conditions
    const double p0 = 101325; //Pa
    const double T0 = 288.150; // K
    const double R = 287.0528; // Nm/kgK
    

    if (Z <= 11000.0){
        const double T0_prime = -0.0065; // K/m
        const double Z0 = 0.0; // m
        T = T + T0 + T0_prime*(Z-Z0); // K
        p = p0*pow(((T0 +T0_prime*(Z-Z0))/T0),(-G0/(R*T0_prime)));
        rho = p/(R*T);
        a = sqrt(1.4*R*T);
        std::cout << "check 1" << std::endl;
    } else if (Z > 11000.0 && Z <= 20000.0){
        // integrate accross the first altitude range
        const double T1 = 216.650;
        const double p1 = p0*pow(((T0 +-0.0065*(11000.0))/T0),(-G0/(R*-0.0065)));
        std::cout << p1 << std::endl;
        // add the portion of the second altitude range
        const double T1_prime = 0.0; // K/m
        const double Z1 = 11000.0; // m
        T = T + T1 + T1_prime*(Z-Z1); // K
        // use the Ti prime = 0 equation
        p = p1*exp(-(G0*(Z-Z1))/(R*T1)); // Pa
        rho = p/(R*T);
        a = sqrt(1.4*R*T);
    } else if (32000.0 >= Z && Z > 20000.0){
        // find T1 and p1
        const double T1 = 216.650;
        const double p1 = p0*pow(((T0 + -0.0065*(11000.0))/T0),(-G0/(R*-0.0065)));
        
        // find T2 and p2
        const double T2 = 216.650;
        const double p2 = p1*exp((-G0*(20000.0-11000.0))/(R*T1));
        
        // add the portion of the third altitude range
        const double T2_prime = 0.001; // K/m
        const double Z2 = 20000.0; // m
        T = T + T2 + T2_prime*(Z-Z2); // K
        // use the Ti prime /= 0 equation
        p = p2*pow(((T2 +T2_prime*(Z-Z2))/T2),(-G0/(R*T2_prime)));
        rho = p/(R*T);
        a = sqrt(1.4*R*T);
    } else if (47000.0 >= Z && Z > 32000.0){
        // find T1 and p1
        const double T1 = 216.650;
        const double p1 = p0*pow(((T0 + -0.0065*(11000.0))/T0),(-G0/(R*-0.0065)));
        
        // find T2 and p2
        const double T2 = 216.650;
        const double p2 = p1*exp((-G0*(20000.0-11000.0))/(R*T1));
        
        // find T3 and p3
        const double T3 = 228.650;
        const double p3 = p2*pow(((T2 + 0.001*(32000.0-20000.0))/T2),(-G0/(R*0.001)));

        // add the portion of the third altitude range
        const double T3_prime = 0.0028; // K/m
        const double Z3 = 32000.0; // m
        T = T + T3 + T3_prime*(Z-Z3); // K
        // use the Ti prime /= 0 equation
        p = p3*pow(((T3 + T3_prime*(Z-Z3))/T3),(-G0/(R*T3_prime)));
        rho = p/(R*T);
        a = sqrt(1.4*R*T);
    } else if (52000.0 >= Z && Z > 47000.0){
        // find T1 and p1
        const double T1 = 216.650;
        const double p1 = p0*pow(((T0 + -0.0065*(11000.0))/T0),(-G0/(R*-0.0065)));
        
        // find T2 and p2
        const double T2 = 216.650;
        const double p2 = p1*exp((-G0*(20000.0-11000.0))/(R*T1));
        
        // find T3 and p3
        const double T3 = 228.650;
        const double p3 = p2*pow(((T2 + 0.001*(32000.0-20000.0))/T2),(-G0/(R*0.001)));

        // find T4 and p4
        const double T4 = 270.650;
        const double p4 = p3*pow(((T3 + 0.0028*(47000.0-32000.0))/T3),(-G0/(R*0.0028)));

        // add the portion of the third altitude range
        const double T4_prime = 0.0; // K/m
        const double Z4 = 47000.0; // m
        T = T + T4 + T4_prime*(Z-Z4); // K
        // use the Ti prime /= 0 equation
        p = p4*exp((-G0*(Z-Z4))/(R*T4));
        rho = p/(R*T);
        a = sqrt(1.4*R*T);
    } else if (61000.0 >= Z && Z > 52000.0){
        // find T1 and p1
        const double T1 = 216.650;
        const double p1 = p0*pow(((T0 + -0.0065*(11000.0))/T0),(-G0/(R*-0.0065)));
        
        // find T2 and p2
        const double T2 = 216.650;
        const double p2 = p1*exp((-G0*(20000.0-11000.0))/(R*T1));
        
        // find T3 and p3
        const double T3 = 228.650;
        const double p3 = p2*pow(((T2 + 0.001*(32000.0-20000.0))/T2),(-G0/(R*0.001)));

        // find T4 and p4
        const double T4 = 270.650;
        const double p4 = p3*pow(((T3 + 0.0028*(47000.0-32000.0))/T3),(-G0/(R*0.0028)));

        // find T5 and p5
        const double T5 = 270.650;
        const double p5 = p4*exp((-G0*(52000.0-47000.0))/(R*T4));

        // add the portion of the third altitude range
        const double T5_prime = -0.002; // K/m
        const double Z5 = 52000.0; // m
        T = T + T5 + T5_prime*(Z-Z5); // K
        // use the Ti prime /= 0 equation
        p = p5*pow(((T5 + T5_prime*(Z-Z5))/T5),(-G0/(R*T5_prime)));
        rho = p/(R*T);
        a = sqrt(1.4*R*T);
    } else if (79000.0 >= Z && Z > 61000.0){
        // find T1 and p1
        const double T1 = 216.650;
        const double p1 = p0*pow(((T0 + -0.0065*(11000.0))/T0),(-G0/(R*-0.0065)));
        
        // find T2 and p2
        const double T2 = 216.650;
        const double p2 = p1*exp((-G0*(20000.0-11000.0))/(R*T1));
        
        // find T3 and p3
        const double T3 = 228.650;
        const double p3 = p2*pow(((T2 + 0.001*(32000.0-20000.0))/T2),(-G0/(R*0.001)));

        // find T4 and p4
        const double T4 = 270.650;
        const double p4 = p3*pow(((T3 + 0.0028*(47000.0-32000.0))/T3),(-G0/(R*0.0028)));

        // find T5 and p5
        const double T5 = 270.650;
        const double p5 = p4*exp((-G0*(52000.0-47000.0))/(R*T4));

        // find T6 and p6
        const double T6 = 252.650;
        const double p6 = p5*pow(((T5 + -0.002*(61000.0-52000.0))/T5),(-G0/(R*-0.002)));

        // add the portion of the third altitude range
        const double T6_prime = -0.004; // K/m
        const double Z6 = 61000.0; // m
        T = T + T6 + T6_prime*(Z-Z6); // K
        // use the Ti prime /= 0 equation
        p = p6*pow(((T6 + T6_prime*(Z-Z6))/T6),(-G0/(R*T6_prime)));
        rho = p/(R*T);
        a = sqrt(1.4*R*T);
    }


    ans[0] = Z;
    ans[1] = T;
    ans[2] = p;
    ans[3] = rho;
    ans[4] = a;
    return ans;
}