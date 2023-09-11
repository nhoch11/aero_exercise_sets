#include "arrow.h"


arrow::arrow(string filename)
{
    std::ifstream f(filename);
    json data = json::parse(f);
    
    // read in json data and assign to existing class attributes
    m_time_step = data["simulation"]["time_step[s]"];
    m_ref_area = data["reference"]["area[ft^2]"];
    m_ref_length = data["reference"]["length[ft]"];
    m_V = data["initial"]["airspeed[ft/s]"];
    m_z = data["initial"]["altitude[ft]"];
    m_theta = data["initial"]["elevation_angle[deg]"];
    m_weight = data["mass"]["weight[lbf]"];
    m_Iyy = data["mass"]["Iyy[slug*ft^2]"];
    m_CL_a = data["aerodynamics"]["CL,a"];
    m_CD0 = data["aerodynamics"]["CD2"];
    m_Cm_a = data["aerodynamics"]["Cm,a"];
    m_Cm_q = data["aerodynamics"]["Cm,q"];

    //y0 = u     dy0 = u_dot
    //y1 = w     dy1 = w_dot
    //y2 = q     dy2 = q_dot
    //y3 = xf    dy3 = xf_dot
    //y4 = zf    dy4 = zf_dot
    //y5 = theta    dy5 = theta_dot
    m_mass = m_weight/gravity_english(m_Z);

    // initialize state ( other indices are zero)
    m_y[0] = m_V;
    m_y[4] = -m_z;

    // dy is all zeros to start with

    
}

void aerodynamics_2_2(double* y, double* FM)
{
    
    // calculate alpha, V
    m_alpha  = atan2(m_y[1], m_y[0]);
    m_V = sqrt(pow(m_y[0],2) + pow(m_Y[1],2));

    // calculate CL, CD, and Cm
    m_CL = m_CL_a*m_alpha;
    m_CD = m_CD0 + m_CD2*pow(m_CL,2);
    m_Cm = m_Cm_a*m_alpha + m_Cm_q*m_ref_length*m_y[2]/m_V;

    // get rho
    get_atmospheric_properties_english(-m_y[4], Atmosphere& atm);
    m_rho = atm.density

    // update F_xb, F_zb, and M_yb
    m_F_xb = -0.5*m_rho*pow(m_V,2)*m_ref_area*m_CD;
    m_F_zb = -0.5*m_rho*pow(m_V,2)*m_ref_area*m_CL;
    m_M_yb = 0.5*m_rho*pow(m_V,2)*m_ref_area*m_ref_length*m_Cm;
}
double* arrow_EOM(double t, double* m_y, double m_mass, double* FM)
{
    m_y[4] = m_mass*gravity_english(m_y[4]);

    // updatae aerodynamic data(call aerodynamics_2_2)
    aerodynamics_2_2(m_y, FM);

    m_g = gravity_english(-m_y[4])

    m_dy[0] = m_g*m_F_xb/(m_mass*m_g) - m_g*sin(m_y[5]) - m_y[2]*m_y[1];
    m_dy[1] = m_g*m_F_zb/(m_mass*m_g) - m_g*cos(m_y[5]) - m_y[2]*m_y[0];
    m_dy[2] = m_M_yb/m_Iyy;
    m_dy[3] = m_y[0]*cos(m_y[5]) + m_y[1]*sin(m_y[5]);
    m_dy[4] = -m_y[0]*sin(m_y[5]) + m_y[1]*cos(m_y[5]);
    m_dy[5] = m_y[2];

    return m_dy;

    
}