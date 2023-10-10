#include "aircraft.h"

// update aerodynamic coeifficint calcs in aerodynamics
// then update forces in func function


aircraft::aircraft(string filename)
{
    std::ifstream f(filename);
    json data = json::parse(f);
    
    // read in json data and assign to existing class attributes
    
    // simulation
    m_time_step = data["simulation"]["time_step[sec]"];
    m_total_time = data["simulation"]["total_time[sec]"];

    // aircraft
    m_wing_area = data["aircraft"]["wing_area[ft^2]"];
    m_wing_span = data["aircraft"]["wing_span[ft]"];
    m_weight = data["aircraft"]["weight[lbf]"];
    m_Ixx = data["aircraft"]["Ixx[slug-ft^2]"];
    m_Iyy = data["aircraft"]["Iyy[slug-ft^2]"];
    m_Izz = data["aircraft"]["Izz[slug-ft^2]"];
    m_Ixy = data["aircraft"]["Ixy[slug-ft^2]"];
    m_Ixz = data["aircraft"]["Ixz[slug-ft^2]"];
    m_Iyz = data["aircraft"]["Iyz[slug-ft^2]"];
    m_hx = data["aircraft"]["hx[slug-ft^2/s]"];
    m_hy = data["aircraft"]["hy[slug-ft^2/s]"];
    m_hz = data["aircraft"]["hz[slug-ft^2/s]"];
    m_cg_shift = data["aircraft"]["CG-shift[ft]"];
    // put thrust stuff here
    m_T0 = data["aircraft"]["thrust"]["T0[lbf]"];
    m_T1 = data["aircraft"]["thrust"]["T1[lbf-s/ft]"];
    m_T2 = data["aircraft"]["thrust"]["T2[lbf-s^2/ft^2]"];
    m_a = data["aircraft"]["thrust"]["a"];

    // initial state
    m_V = data["initial"]["airspeed[ft/s]"];
    m_altitude = data["initial"]["altitude[ft]"];
    m_elv = data["initial"]["elevation_angle[deg]"]*180.0/pi;
    m_bank = data["initial"]["bank_angle[deg]"]*180.0/pi;
    m_alpha = data["initial"]["alpha[deg]"]*180.0/pi;
    m_beta = data["initial"]["beta[deg]"]*180.0/pi;
    m_p = data["initial"]["p[deg/s]"];
    m_q = data["initial"]["q[deg/s]"];
    m_r = data["initial"]["r[deg/s]"];
    m_heading = data["initial"]["heading_angle[deg]"]*180.0/pi;
    m_aileron = data["initial"]["aileron[deg]"]*180.0/pi;
    m_elevator = data["initial"]["elevator[deg]"]*180.0/pi;
    m_rudder = data["initial"]["rudder[deg]"]*180.0/pi;
    m_throttle = data["initial"]["throttle"];

    // aerodynamics
    m_CL0 = data["aerodynamics"]["CL0"];
    m_CL_a = data["aerodynamics"]["CL,a"];
    m_CL_qbar = data["aerodynamics"]["CL,qbar"];
    m_CL_de = data["aerodynamics"]["CL,de"];
    m_CS_b = data["aerodynamics"]["CS,b"];
    m_CS_pbar = data["aerodynamics"]["CS,pbar"];
    m_CS_rbar = data["aerodynamics"]["CS,rbar"];
    m_CS_da = data["aerodynamics"]["CS,da"];
    m_CS_dr = data["aerodynamics"]["CS,dr"];
    m_CDL0 = data["aerodynamics"]["CDL0"];
    m_CD_L = data["aerodynamics"]["CD,L"];
    m_CD_L2 = data["aerodynamics"]["CD,L2"];
    m_CD_S2 = data["aerodynamics"]["CD,S2"];
    m_CD_qbar = data["aerodynamics"]["CD,qbar"];
    m_CD_Lqbar = data["aerodynamics"]["CD,Lqbar"];
    m_CD_de = data["aerodynamics"]["CD,de"];
    m_CD_Lde = data["aerodynamics"]["CD,Lde"];
    m_CD_de2 = data["aerodynamics"]["CD,de2"];
    m_Cl_b = data["aerodynamics"]["Cl,b"];
    m_Cl_pbar = data["aerodynamics"]["Cl,pbar"];
    m_Cl_rbar = data["aerodynamics"]["Cl,rbar"];
    m_Cl_Lrbar = data["aerodynamics"]["Cl,Lrbar"];
    m_Cl_da = data["aerodynamics"]["Cl,da"];
    m_Cl_dr = data["aerodynamics"]["Cl,dr"];
    m_Cm0 = data["aerodynamics"]["Cm0"];
    m_Cm_a = data["aerodynamics"]["Cm,a"];
    m_Cm_qbar = data["aerodynamics"]["Cm,qbar"];
    m_Cm_de = data["aerodynamics"]["Cm,de"];
    m_Cn_b = data["aerodynamics"]["Cn,b"];
    m_Cn_pbar = data["aerodynamics"]["Cn,pbar"];
    m_Cn_Lpbar = data["aerodynamics"]["Cn,Lpbar"];
    m_Cn_rbar = data["aerodynamics"]["Cn,rbar"];
    m_Cn_da = data["aerodynamics"]["Cn,da"];
    m_Cn_Lda = data["aerodynamics"]["Cn,Lda"];
    m_Cn_dr = data["aerodynamics"]["Cn,dr"];
    
}


void aircraft::aerodynamics_aircraft(double* y, double* ans)
{
    // make dummy variables for readibility
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
    //double CL =  m_CLa * alpha;
    //double CS =  m_CLa * beta;
    //double CD =  m_CD0 + (m_CD2 * pow(CL,2));
    double Cl =  m_Cl0 + ((m_Clp * m_ref_length * p) / V);
    double Cm =  (m_Cma * alpha) + ((m_Cmq * m_ref_length * q) / V);
    double Cn = -(m_Cma * beta) +  ((m_Cmq * m_ref_length * r) / V);

    // get rho, mu, Re
    get_atmospheric_properties_english(-zf, m_atm);
    double rho = m_atm.density;
    double mu  = m_atm.dynamic_viscosity;
    double Re  = rho*V*m_ref_length/mu;

    // get drag of a sphere
    double CD = get_sphere_CD(Re);
    // update forces and moments
    ans[0] =  -0.5*rho*pow(V,2)*m_ref_area*CD*cos(alpha)*cos(beta); // F_xb
    ans[1] =  -0.5*rho*pow(V,2)*m_ref_area*CD*sin(beta); // F_yb
    ans[2] =  -0.5*rho*pow(V,2)*m_ref_area*CD*sin(alpha)*cos(beta); // F_zb
    ans[3] =  0.0; // 0.5*rho*pow(V,2)*m_ref_area*m_ref_length*Cl; // M_xb
    ans[4] =  0.0; // 0.5*rho*pow(V,2)*m_ref_area*m_ref_length*Cm; // M_yb
    ans[5] =  0.0; // 0.5*rho*pow(V,2)*m_ref_area*m_ref_length*Cn; // M_zb
}

void aircraft::aircraft_rk4_func(double t, double* y, double* ans)
{

    double* FM = new double[6];
    // update aerodynamic data
    aerodynamics_aircraft(y, FM);
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
    double e0 = y[9];
    double ex = y[10];
    double ey = y[11];
    double ez = y[12];

    double g = gravity_english(-zf);

    ans[0]  = (g*Fxb/m_weight) + (g*2.0*(ex*ez - ey*e0)) + (r*v) - (q*w); // udot
    ans[1]  = (g*Fyb/m_weight) + (g*2.0*(ey*ez + ex*e0)) + (p*w) - (r*u); // vdot
    ans[2]  = (g*Fzb/m_weight) + (g*(ez*ez + e0*e0 - ex*ex - ey*ey)) + (q*u) - (p*v); // wdot
    ans[3]  = (Mxb/m_Ixx); // pdot
    ans[4]  = (Myb + (m_Iyy - m_Ixx)*p*r)/m_Iyy; // qdot
    ans[5]  = (Mzb + (m_Ixx - m_Iyy)*p*q)/m_Iyy; // rdot
    
    // build quat vectors for first quat mult
    double* quatA    = new double[4];
    double* quatB    = new double[4];
    double* quat_AB  = new double[4];
    double* quat_xyz = new double[4];
    
    quatA[0] = 0.0;
    quatA[1] = u;
    quatA[2] = v;
    quatA[3] = w;

    quatB[0] =  e0;
    quatB[1] = -ex;
    quatB[2] = -ey;
    quatB[3] = -ez;

    quat_mult(quatA, quatB, quat_AB);

    quat_mult(&y[9], quat_AB, quat_xyz);

    // take the last 3 components of quat_xyz as x, y and z (the first element is 0.0)
    ans[6] = quat_xyz[1];
    ans[7] = quat_xyz[2];
    ans[8] = quat_xyz[3];

    ans[9]  = 0.5*(-ex*p - ey*q - ez*r); // e0 dot
    ans[10] = 0.5*( e0*p - ez*q + ey*r); // ex dot
    ans[11] = 0.5*( ez*p + e0*q - ex*r); // ey dot
    ans[12] = 0.5*(-ey*p + ex*q + e0*r); // ez dot
    
}

void aircraft::aircraft_rk4(double t0, double* y0, double dt, int size, double* ans)
{
    // print first function calls in check file
    FILE* check_file = fopen("check_4_3.txt", "w");
    fprintf(check_file, "  Time[s]              udot                vdot                 wdot                  pdot               qdot                 rdot                 xdot                  ydot                 zdot                   e0dot                exdot               eydot              ezdot\n");
    
    double* dy = new double[size];
    // k1:
    aircraft_rk4_func(t0, y0, dy);
    fprintf(check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11], dy[12]);
    
    // multiply k1 by dt 
    for (int i = 0; i < size; i++)
    {
        m_k1[i] = dt * dy[i];
        // use k1 to update y_temporary
        m_y_temp[i] = y0[i] + 0.5 * m_k1[i];
    }

    // k2:  
    aircraft_rk4_func(t0 + (0.5 * dt), m_y_temp, dy);
    fprintf(check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11], dy[12]);
    
    // multiply k2 by dt and then update y_temp
    for (int i = 0; i < size; i++)
    {
        m_k2[i] = dt * dy[i];
        m_y_temp[i] = y0[i] + (0.5 * m_k2[i]);
    }
    
    // k3:
    aircraft_rk4_func(t0 + (0.5 * dt), m_y_temp, dy);
    fprintf(check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11], dy[12]);//fprintf(check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11]);
    
    // multiply k2 by dt and then update y_temp
    for (int i = 0; i < size; i++)
    {  
        m_k3[i] = dt * dy[i];
        m_y_temp[i] = y0[i] + m_k3[i];
    }
    
    // k4:
    aircraft_rk4_func(t0 + dt, m_y_temp, dy);
    fprintf(check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11], dy[12]);//fprintf(check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11]);
    fclose(check_file);
    
    // multiply k4 by dt   
    for (int i = 0; i < size; i++)
    {
        m_k4[i] = dt * dy[i];

        // find new
        ans[i] = y0[i] + (m_k1[i] + 2.0 * m_k2[i] + 2.0 * m_k3[i] + m_k4[i]) / 6.0  ;
    }
    
}


void aircraft::init_sim()
{
    Atmosphere atm;
    
    // print results in a file
    FILE* en_file = fopen("aircraft.txt", "w");
    fprintf(en_file, "  Time[s]              u[ft/s]            v[ft/s]                w[ft/s]              p[rad/s]             q[rad/s]             r[rad/s]             x[ft]                y[ft]                z[ft]                e0                   ex                   ey                   ez\n");
    
    double t0 = 0.0;
    int size = 13;
    double* y0 = new double[size];
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

    // convert phi, theta, and psi to a quat
    double* quat = new double[4];
    euler_to_quat(&y0[9], quat);

    // normalize the quat
    quat_norm(quat);

    // store the quat in y0
    y0[9] = quat[0];
    y0[10] = quat[1];
    y0[11] = quat[2];
    y0[12] = quat[3];

    double* y = new double[13];
}

void aircraft::exercise_6_1()
{
    
    do {
        aircraft_rk4(t0, y0, m_time_step, size, y);
        
        fprintf(en_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, y0[0], y0[1], y0[2],y0[3],y0[4], y0[5], y0[6], y0[7], y0[8], y0[9], y0[10], y0[11], y0[12]);
        
        // y = y0 ( copy y0 into y)
        array_copy(y, y0, size);

        // re-normalize the quat
        quat_norm(&y[9]);
        
        // add time step
        t0 += m_time_step;
    } while (-y0[8] > 0.0);
    //while (t0<0.005);
    
    fclose(en_file);
    
}

