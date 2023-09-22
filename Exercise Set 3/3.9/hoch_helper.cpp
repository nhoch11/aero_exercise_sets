#include "hoch_helper.h"

void get_atmospheric_properties_si(double altitude, Atmosphere& atm)
{
    // declare constants and initialize variables
    const double RE = 6356766.0; // meters
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
   
   // first layer
    if (Z <= 11000.0){
        T = T_array[0] + T_prime[0]*(Z-Z_const[0]); // K
        p = p_array[0]*pow(((T_array[0] +T_prime[0]*(Z-Z_const[0]))/T_array[0]),(-G0/(R*T_prime[0])));
    } 
    
    // second layer
    else if (Z > 11000.0 && Z <= 20000.0){
       // add the portion of the second altitude range to the previously integrated first layer
        T = T_array[1] + T_prime[1]*(Z-Z_const[1]); // K
        // use the Ti prime = 0 equation
        p = p_array[1]*exp(-(G0*(Z-Z_const[1]))/(R*T_array[1])); // Pa
    
    } 
    
    // third layer
    else if (32000.0 >= Z && Z > 20000.0){
        // add the portion of the third altitude range to the previously integrated 2 layers
        T = T_array[2] + T_prime[2]*(Z-Z_const[2]); // K
        // use the Ti prime /= 0 equation
        p = p_array[2]*pow(((T_array[2] +T_prime[2]*(Z-Z_const[2]))/T_array[2]),(-G0/(R*T_prime[2])));
        
    } 
    
    //fourth layer
    else if (47000.0 >= Z && Z > 32000.0){
        // add the portion of the fourth altitude to the previously integrated 3 layers
        T = T_array[3] + T_prime[3]*(Z-Z_const[3]); // K
        // use the Ti prime /= 0 equation
        p = p_array[3]*pow(((T_array[3] + T_prime[3]*(Z-Z_const[3]))/T_array[3]),(-G0/(R*T_prime[3])));
        
    } 
    
    // fifth layer
    else if (52000.0 >= Z && Z > 47000.0){
        // add the portion of the fifth altitude range to the previously integrated 4 layers
        T = T_array[4] + T_prime[4]*(Z-Z_const[4]); // K
        // use the Ti prime /= 0 equation
        p = p_array[4]*exp((-G0*(Z-Z_const[4]))/(R*T_array[4]));
        
    } 
    
    // sixth layer
    else if (61000.0 >= Z && Z > 52000.0){
        // add the portion of the sixth altitude range to the prevously integrated 5 layers
        T = T_array[5] + T_prime[5]*(Z-Z_const[5]); // K
        // use the Ti prime /= 0 equation
        p = p_array[5]*pow(((T_array[5]  + T_prime[5]*(Z-Z_const[5]))/T_array[5] ),(-G0/(R*T_prime[5])));
        
    } 
    
    // seventh layer
    else if (79000.0 >= Z && Z > 61000.0){
        // add the portion of the seventh altitude range to the previously integrated 6 layers
        T = T_array[6] + T_prime[6]*(Z-Z_const[6]); // K
        // use the Ti prime /= 0 equation
        p = p_array[6]*pow(((T_array[6] + T_prime[6]*(Z-Z_const[6]))/ T_array[6]),(-G0/(R*T_prime[6])));
        
    }
    rho = p/(R*T);
    a = sqrt(1.4*R*T);

    atm.geopotential_altitude = Z;
    atm.temperature = T;
    atm.pressure = p;
    atm.density = rho;
    atm.speed_of_sound = a;
}

void get_atmospheric_properties_english(double altitude, Atmosphere& atm)
{
    // convert input of feet to meters
    altitude = altitude*0.3048;

    // call get atmosphere si 
    get_atmospheric_properties_si(altitude, atm);

    // convert results from si to english
    atm.geopotential_altitude = atm.geopotential_altitude/0.3048; // meters to ft
    atm.temperature = atm.temperature*1.8;  // Kelvin to Rankine
    atm.pressure = atm.pressure*0.020885434304801722; // Pa to lbf/ft^2
    atm.density = atm.density*0.00194032032363104; // kg/m^3 to slugs/ ft^3
    atm.speed_of_sound = atm.speed_of_sound/0.3048; // m/s to ft/s
    
}

void rk4_array(double t0, double* y0, double dt, int size, double* ans)
{
    // declare variables for k's
    double* k1 = new double[size];
    double* k2 = new double[size];
    double* k3 = new double[size];
    double* k4 = new double[size];

    double* y_temp = new double[size];
    
    // k1:
    k1 = f_array(t0, y0);

    // multiply k1 by dt 
    for (int i = 0; i < size; i++)
    {
        k1[i] = dt * k1[i];
        // use k1 to update y_temporary
        y_temp[i] = y0[i] + 0.5*k1[i];
    }

    // k2:  
    k2 = f_array(t0 + 0.5*dt, y_temp);
    
    // multiply k2 by dt and then update y_temp
    for (int i = 0; i < size; i++)
    {
        k2[i] = dt * k2[i];
        y_temp[i] = y0[i] + 0.5*k2[i];
    }
    
    // k3:
    k3 = f_array(t0 + 0.5*dt, y_temp);
    
    // multiply k2 by dt and then update y_temp
    for (int i = 0; i < size; i++)
    {  
        k3[i] = dt * k3[i];
        y_temp[i] = y0[i] + k3[i];
    }
    
    // k4:
    k4 = f_array(t0 + dt, y_temp);

    // multiply k4 by dt   
    for (int i = 0; i < size; i++)
    {
        k4[i] = dt * k4[i];

        // find dxdt, then dzdt
        ans[i] = y0[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0  ;
    }
    
    
}

void array_copy(double* A, double* B, int size)
{
    // B = A ( copy A into B)
    for (int i=0; i < size; i++){
        B[i] = A[i];
    }
}

void quat_mult(double* A, double* B, double* ans)
{
    double a0 = A[0];
    double ax = A[1];
    double ay = A[2];
    double az = A[3];
    double b0 = B[0];
    double bx = B[1];
    double by = B[2];
    double bz = B[3];
    // equation 11.6.3
    ans[0] = a0*b0 - ax*bx - ay*by - az*bz;
    ans[1] = a0*bx + ax*b0 + ay*bz - az*by;
    ans[2] = a0*by - ax*bz + ay*b0 + az*bx;
    ans[3] = a0*bz + ax*by - ay*bx + az*b0;
}

void quat_norm(double* quat)
{   // equation 11.6.5
    double magnitude = sqrt(pow(quat[0], 2) + pow(quat[1], 2) + pow(quat[2], 2) + pow(quat[3], 2));
    quat[0] /= magnitude;
    quat[1] /= magnitude;
    quat[2] /= magnitude;
    quat[3] /= magnitude;
}

void euler_to_quat(double* eul, double* quat)
{ // equation 11.7.8
    double cphi2   = cos(eul[0]/2.0);
    double ctheta2 = cos(eul[1]/2.0);
    double cpsi2   = cos(eul[2]/2.0);
    double sphi2   = sin(eul[0]/2.0);
    double stheta2 = sin(eul[1]/2.0);
    double spsi2   = sin(eul[2]/2.0);
    quat[0] = cphi2*ctheta2*cpsi2 + sphi2*stheta2*spsi2;
    quat[1] = sphi2*ctheta2*cpsi2 - cphi2*stheta2*spsi2;
    quat[2] = cphi2*stheta2*cpsi2 + sphi2*ctheta2*spsi2;
    quat[3] = cphi2*ctheta2*spsi2 - sphi2*stheta2*cpsi2;
}

void quat_to_euler(double* quat, double* eul)
{ // equation 11.7.11
    double e0 = quat[0];
    double ex = quat[1];
    double ey = quat[2];
    double ez = quat[3];

    if (e0*ey - ex*ez == 0.5)
    {
        eul[0] = 2.0*asin(ex/cos(pi/4.0));
        eul[1] = pi/2.0;
        eul[2] = 0.0;
    }

    if (e0*ey - ex*ez == -0.5)
    {
        eul[0] = 2.0*asin(ex/cos(pi/4.0));
        eul[1] = -pi/2.0;
        eul[2] = 0.0;
    }

    else
    {
        eul[0] = atan2( 2.0*(e0*ex + ey*ez), pow(e0, 2.0) + pow(ez, 2.0) - pow(ex, 2.0) - pow(ey, 2.0));
        eul[1] = asin(2.0*(e0*ey - ex*ez));
        eul[2] = atan2( 2.0*(e0*ez + ex*ey), pow(e0, 2.0) + pow(ex, 2.0) - pow(ey, 2.0) - pow(ez, 2.0));
    }
}

void body_to_fixed(double* vec, double* quat, double* ans)
{   
    double vxb = vec[0];
    double vyb = vec[1];
    double vzb = vec[2];
    double e0 = quat[0];
    double ex = quat[1];
    double ey = quat[2];
    double ez = quat[3];
    // equation 11.5.6

    ans[0] = (pow(ex,2.0) + pow(e0,2.0) - pow(ey,2.0) - pow(ez,2.0))*vxb + 2.0*(ex*ey - ez*e0)*vyb + 2.0*(ex*ez + ey*e0)*vzb; 
    ans[1] = 2.0*(ex*ey + ez*e0)*vxb + (pow(ey,2.0) + pow(e0,2.0) - pow(ex,2.0) - pow(ez,2.0))*vyb + 2.0*(ey*ez - ex*e0)*vzb;
    ans[2] = 2.0*(ex*ez - ey*e0)*vxb + 2.0*(ey*ez + ex*e0)*vyb + (pow(ez,2.0) + pow(e0,2.0) - pow(ex,2.0) - pow(ey,2.0))*vzb;

}

void fixed_to_body(double* vec, double* quat, double* ans)
{
    double vxb = vec[0];
    double vyb = vec[1];
    double vzb = vec[2];
    double e0 = quat[0];
    double ex = quat[1];
    double ey = quat[2];
    double ez = quat[3];
    // equation 11.5.5

    ans[0] = (pow(ex,2.0) + pow(e0,2.0) - pow(ey,2.0) - pow(ez,2.0))*vxb + 2.0*(ex*ey + ez*e0)*vyb + 2.0*(ex*ez - ey*e0)*vzb; 
    ans[1] = 2.0*(ex*ey - ez*e0)*vxb + (pow(ey,2.0) + pow(e0,2.0) - pow(ex,2.0) - pow(ez,2.0))*vyb + 2.0*(ey*ez + ex*e0)*vzb;
    ans[2] = 2.0*(ex*ez + ey*e0)*vxb + 2.0*(ey*ez - ex*e0)*vyb + (pow(ez,2.0) + pow(e0,2.0) - pow(ex,2.0) - pow(ey,2.0))*vzb;

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

double gravity_si(double H)
{
    // H given in meters
    double g = 9.806645*pow((6356766.0/ (6356766.0 + H)),2); // m/s^2
    return g;
}

double gravity_english(double H)
{
    // convert H feet to H meters
    H =  H*0.3048;

    // call gravity si and convert to ft/s^2

    double g = gravity_si(H)*3.28083989501312 ;

    return g;
}

