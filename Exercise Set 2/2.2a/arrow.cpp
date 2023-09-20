#include "arrow.h"


arrow::arrow(string filename)
{
    std::ifstream f(filename);
    json data = json::parse(f);
    
    // read in json data and assign to existing class attributes
    m_time_step = data["simulation"]["time_step[s]"];
    m_ref_area = data["reference"]["area[ft^2]"];
    m_ref_length = data["reference"]["length[ft]"];
    m_init_V = data["initial"]["airspeed[ft/s]"];
    m_init_altitude = data["initial"]["altitude[ft]"];
    double theta_deg = data["initial"]["elevation_angle[deg]"];
    m_init_theta = (theta_deg*pi)/180.0;
    m_weight = data["mass"]["weight[lbf]"];
    m_Iyy = data["mass"]["Iyy[slug*ft^2]"];
    m_CLa = data["aerodynamics"]["CL,a"];
    m_CD0 = data["aerodynamics"]["CD0"];
    m_CD2 = data["aerodynamics"]["CD2"];
    m_Cma = data["aerodynamics"]["Cm,a"];
    m_Cmq = data["aerodynamics"]["Cm,q"];

    

    //y0 = u     dy0 = u_dot
    //y1 = w     dy1 = w_dot
    //y2 = q     dy2 = q_dot
    //y3 = xf    dy3 = xf_dot
    //y4 = zf    dy4 = zf_dot
    //y5 = theta    dy5 = theta_dot

    // initialize state ( other indices are zero)
    // m_y[0] = m_V;
    // m_y[4] = -m_init_altitude;

    // dy is all zeros to start with

    //get_atmospheric_properties_english(-m_init_altitude, m_atm);
    
    
}

void arrow::aerodynamics_2_2(double* y, double* ans)
{
    //double* ans = new double[3];
    // calculate alpha, V
    double alpha  = atan2(y[1], y[0]);
    double V = sqrt(pow(y[0],2) + pow(y[1],2));

    // calculate CL, CD, and Cm
    double CL = m_CLa*alpha;
    double CD = m_CD0 + m_CD2*pow(CL,2);
    double Cm = m_Cma*alpha + (m_Cmq*m_ref_length*y[2])/V;

    // get rho
    get_atmospheric_properties_english(-y[4], m_atm);
    double rho = m_atm.density;

    // update F_xb [0], F_zb [1], and M_yb [2]
    ans[0] = -0.5*rho*pow(V,2)*m_ref_area*CD;
    ans[1] = -0.5*rho*pow(V,2)*m_ref_area*CL;
    ans[2] = 0.5*rho*pow(V,2)*m_ref_area*m_ref_length*Cm;
}

void arrow::arrow_rk4_func_2_2(double t, double* y, double* ans)
{
    
    double* FM = new double[3];
    // updatae aerodynamic data(call aerodynamics_2_2)
    aerodynamics_2_2(y, FM);

    double g = gravity_english(-y[4]);

    ans[0] = (g*FM[0])/(m_weight) - g*sin(y[5]) - y[2]*y[1];
    ans[1] = (g*FM[1])/(m_weight) + g*cos(y[5]) + y[2]*y[0];
    ans[2] = FM[2]/m_Iyy;
    ans[3] = y[0]*cos(y[5]) + y[1]*sin(y[5]);
    ans[4] = -y[0]*sin(y[5]) + y[1]*cos(y[5]);
    ans[5] = y[2];
    
}

void arrow::arrow_rk4_2_2(double t0, double* y0, double dt, int size, double* ans)
{
    double* dy = new double[6];
    // k1:
    arrow_rk4_func_2_2(t0, y0, dy);
    //m_k1 = dy;

    // multiply k1 by dt 
    for (int i = 0; i < size; i++)
    {
        m_k1[i] = dt * dy[i];
        // use k1 to update y_temporary
        m_y_temp[i] = y0[i] + 0.5*m_k1[i];
    }

    // k2:  
    arrow_rk4_func_2_2(t0 + 0.5*dt, m_y_temp, dy);
    //m_k2 = dy;
    
    // multiply k2 by dt and then update y_temp
    for (int i = 0; i < size; i++)
    {
        m_k2[i] = dt * dy[i];
        m_y_temp[i] = y0[i] + 0.5*m_k2[i];
    }
    
    // k3:
    arrow_rk4_func_2_2(t0 + 0.5*dt, m_y_temp, dy);
    //m_k3 = dy;
    // multiply k2 by dt and then update y_temp
    for (int i = 0; i < size; i++)
    {  
        m_k3[i] = dt * dy[i];
        m_y_temp[i] = y0[i] + m_k3[i];
    }
    
    // k4:
    arrow_rk4_func_2_2(t0 + dt, m_y_temp, dy);
    //m_k4 = dy;

    // multiply k4 by dt   
    for (int i = 0; i < size; i++)
    {
        m_k4[i] = dt * dy[i];

        // find new
        ans[i] = y0[i] + (m_k1[i] + 2.0*m_k2[i] + 2.0*m_k3[i] + m_k4[i])/6.0  ;
    }
    
}

void arrow::exercise_2_2()
{
    Atmosphere atm;
    
    FILE* en_file = fopen("arrow_2_2.txt", "w");

    fprintf(en_file, "  Time[s]              u[ft/s]             w[ft/s]               q[rad/s]             x[ft]                z[ft]                theta[rad]\n");
    double t0 = 0.0;
    double* y0 = new double[6];
    y0[0] = m_init_V;
    y0[1] = 0.0;
    y0[4] = -m_init_altitude;
    y0[5] = m_init_theta;
    int size = 6;
    double* y = new double[6];
    //for (t0 = 0.0; t0 < 1.35; t0 += m_time_step) {
    do {
        arrow_rk4_2_2(t0, y0, m_time_step, size, y);
        fprintf(en_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, y0[0], y0[1], y0[2],y0[3],y0[4], y0[5]);
        y0[0] = y[0];
        y0[1] = y[1];
        y0[2] = y[2];
        y0[3] = y[3];
        y0[4] = y[4];
        y0[5] = y[5];
        t0 += m_time_step;
    } while (-y0[4] > 0.0);
    fclose(en_file);
    
}