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
    m_Ixx = data["mass"]["Ixx[slug*ft^2]"];
    m_Iyy = data["mass"]["Iyy[slug*ft^2]"];
    m_CLa = data["aerodynamics"]["CL,a"];
    m_CD0 = data["aerodynamics"]["CD0"];
    m_CD2 = data["aerodynamics"]["CD2"];
    m_Cma = data["aerodynamics"]["Cm,a"];
    m_Cmq = data["aerodynamics"]["Cm,q"];
    m_Clp = data["aerodynamics"]["Cl,p"];
    m_Cl0 = data["aerodynamics"]["Cl0"];

    

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
    double* dy = new double[size];
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

void arrow::aerodynamics_2_4(double* y, double* ans)
{
    double u  = y[0];
    double v  = y[1];
    double w  = y[2];
    double p  = y[3]; 
    double q  = y[4];
    double r  = y[5];
    double xf = y[6]; 
    double yf = y[7];
    double zf = y[8];
    double phi = y[9];
    double theta = y[10];
    double psi = y[11];
 
    // calculate alpha, V
    double alpha  = atan2(w, u);
    double beta = atan2(v,u);
    double V = sqrt(pow(u,2) + pow(v,2) + pow(w,2));

    // calculate CL, CD, and Cm
    double CL =  m_CLa * alpha;
    double CS =  m_CLa * beta;
    double CD =  m_CD0 + (m_CD2 * pow(CL,2));
    double Cl =  m_Cl0 + ((m_Clp * m_ref_length * p) / V);
    double Cm =  (m_Cma * alpha) + ((m_Cmq * m_ref_length * q) / V);
    double Cn = -(m_Cma * beta) +  ((m_Cmq * m_ref_length * r) / V);

    // get rho
    get_atmospheric_properties_english(-zf, m_atm);
    double rho = m_atm.density;

    // update forces and moments
    ans[0] = -0.5 * rho * pow(V,2) * m_ref_area * CD; // F_xb
    ans[1] = -0.5 * rho * pow(V,2) * m_ref_area * CS; // F_yb
    ans[2] = -0.5 * rho * pow(V,2) * m_ref_area * CL; // F_zb
    ans[3] =  0.5 * rho * pow(V,2) * m_ref_area * m_ref_length * Cl; // M_xb
    ans[4] =  0.5 * rho * pow(V,2) * m_ref_area * m_ref_length * Cm; // M_yb
    ans[5] =  0.5 * rho * pow(V,2) * m_ref_area * m_ref_length * Cn; // M_zb
}

void arrow::arrow_rk4_func_2_4(double t, double* y, double* ans)
{
    
    
    double* FM = new double[6];
    // updatae aerodynamic data(call aerodynamics_2_2)
    aerodynamics_2_4(y, FM);
    double Fxb = FM[0];
    double Fyb = FM[1];
    double Fzb = FM[2];
    double Mxb = FM[3];
    double Myb = FM[4];
    double Mzb = FM[5];

    // declare variables to keep track of stuff
    double u  = y[0];
    double v  = y[1];
    double w  = y[2];
    double p  = y[3]; 
    double q  = y[4];
    double r  = y[5];
    double xf = y[6]; 
    double yf = y[7];
    double zf = y[8];
    double phi = y[9];
    double theta = y[10];
    double psi = y[11];

    double g = gravity_english(-zf);

    ans[0]  = (g * Fxb / m_weight) - (g * sin(theta)) + (r * v) - (q * w); // udot
    ans[1]  = (g * Fyb / m_weight) + (g * sin(phi) * cos(theta)) + (p * w) - (r * u); // vdot
    ans[2]  = (g * Fzb / m_weight) + (g * cos(phi) * cos(theta)) + (q * u) - (p * v); // wdot
    ans[3]  = (Mxb / m_Ixx); // pdot
    ans[4]  = (Myb + (m_Iyy - m_Ixx) * p * r) / m_Iyy; // qdot
    ans[5]  = (Mzb + (m_Ixx - m_Iyy) * p * q) / m_Iyy; // rdot
    ans[6]  = cos(theta) * cos(psi) * u + (sin(phi) * sin(theta) * cos(psi) - cos(phi) * sin(psi)) * v + (cos(phi) * sin(theta) * cos(psi) + sin(phi) * sin(psi)) * w; // xfdot
    ans[7]  = (cos(theta) * sin(psi) * u) + (((sin(phi) * sin(theta) * sin(psi)) + (cos(phi) * cos(psi))) * v) + (((cos(phi) * sin(theta) * sin(psi)) - (sin(phi) * cos(psi))) * w); // yfdot
    ans[8]  = (-sin(theta) * u) + (sin(phi) * cos(theta) * v) + (cos(phi) * cos(theta) * w); // zfdot
    ans[9]  = p + (sin(phi)*sin(theta)/cos(theta))*q + (cos(phi)*sin(theta)/cos(theta))*r; // phi dot
    ans[10] = cos(phi)*q - sin(phi)*r; // theta dot
    ans[11] = ((sin(phi) / cos(theta)) * q) + ((cos(phi) / cos(theta)) * r); // psi dot
    
}

void arrow::arrow_rk4_2_4(double t0, double* y0, double dt, int size, double* ans)
{
    FILE* check_file = fopen("check.txt", "w");

    fprintf(check_file, "  Time[s]              u[ft/s]             v[ft/s]               w[ft/s]               p[rad/s]             q[rad/s]             r[rad/s]             x[ft]                y[ft]                z[ft]                phi[rad]             theta[rad]             psi[rad]\n");
    
    double* dy = new double[size];
    // k1:
    arrow_rk4_func_2_4(t0, y0, dy);
    //m_k1 = dy;
    fprintf(check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11]);
    
    // multiply k1 by dt 
    for (int i = 0; i < size; i++)
    {
        m_k1[i] = dt * dy[i];
        // use k1 to update y_temporary
        m_y_temp[i] = y0[i] + 0.5 * m_k1[i];
    }

    // k2:  
    arrow_rk4_func_2_4(t0 + (0.5 * dt), m_y_temp, dy);
    //m_k2 = dy;
    fprintf(check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11]);
    
    // multiply k2 by dt and then update y_temp
    for (int i = 0; i < size; i++)
    {
        m_k2[i] = dt * dy[i];
        m_y_temp[i] = y0[i] + (0.5 * m_k2[i]);
    }
    
    // k3:
    arrow_rk4_func_2_4(t0 + (0.5 * dt), m_y_temp, dy);
    //m_k3 = dy;
    fprintf(check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11]);
    
    // multiply k2 by dt and then update y_temp
    for (int i = 0; i < size; i++)
    {  
        m_k3[i] = dt * dy[i];
        m_y_temp[i] = y0[i] + m_k3[i];
    }
    
    // k4:
    arrow_rk4_func_2_4(t0 + dt, m_y_temp, dy);
    //m_k4 = dy;
    fprintf(check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11]);
    
    // multiply k4 by dt   
    for (int i = 0; i < size; i++)
    {
        m_k4[i] = dt * dy[i];

        // find new
        ans[i] = y0[i] + (m_k1[i] + 2.0 * m_k2[i] + 2.0 * m_k3[i] + m_k4[i]) / 6.0  ;
    }
    
    fclose(check_file);
}

void arrow::exercise_2_4()
{
    Atmosphere atm;
    
    FILE* en_file = fopen("arrow_2_4.txt", "w");

    fprintf(en_file, "  Time[s]              u[ft/s]             v[ft/s]               w[ft/s]               p[rad/s]             q[rad/s]             r[rad/s]             x[ft]                y[ft]                z[ft]                phi[rad]             theta[rad]             psi[rad]\n");
    double t0 = 0.0;
    double* y0 = new double[12];
    y0[0]  = m_init_V;         // u
    y0[1]  = 0.0;              // v
    y0[2]  = 0.0;              // w
    y0[3]  = 0.0;              // p
    y0[4]  = 0.0;              // q
    y0[5]  = 0.0;              // r
    y0[6]  = 0.0;              // xf
    y0[7]  = 0.0;              // yf
    y0[8]  = -m_init_altitude; // zf
    y0[9]  = 0.0;              // phi
    y0[10] = m_init_theta;     // theta
    y0[11] = 0.0;              // psi
    
    int size = 12;
    double* y = new double[12];
    
    do {
        arrow_rk4_2_4(t0, y0, m_time_step, size, y);
        fprintf(en_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, y0[0], y0[1], y0[2],y0[3],y0[4], y0[5], y0[6], y0[7], y0[8], y0[9], y0[10], y0[11]);
        //for (int i = 0; i<12; i++)
            //{printf("%20.12e", y0[i]);}
        array_copy(y, y0, size);
        
        t0 += m_time_step;
    } 
    //while (t0<0.005);
    while (-y0[8] > 0.0);
    
    fclose(en_file);
    
}