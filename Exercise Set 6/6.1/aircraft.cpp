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
    //m_cg_shift = data["aircraft"]["CG-shift[ft]"];
    // put thrust stuff here
    m_T0 = data["aircraft"]["thrust"]["T0[lbf]"];
    m_T1 = data["aircraft"]["thrust"]["T1[lbf-s/ft]"];
    m_T2 = data["aircraft"]["thrust"]["T2[lbf-s^2/ft^2]"];
    m_a = data["aircraft"]["thrust"]["a"];

    // initial state
    m_V = data["initial"]["airspeed[ft/s]"];
    m_altitude = data["initial"]["altitude[ft]"];
    m_elv_angle = data["initial"]["elevation_angle[deg]"];
    m_elv_angle = m_elv_angle*pi/180.0;
    m_bank = data["initial"]["bank_angle[deg]"];
    m_bank = m_bank*pi/180.0;
    m_alpha = data["initial"]["alpha[deg]"];
    m_alpha = m_alpha*pi/180.0;
    m_beta = data["initial"]["beta[deg]"];
    m_beta = m_beta*pi/180.0;
    m_p = data["initial"]["p[deg/s]"];
    m_p = m_p*pi/180.0;
    m_q = data["initial"]["q[deg/s]"];
    m_q = m_q*pi/180.0;
    m_r = data["initial"]["r[deg/s]"];
    m_r = m_r*pi/180.0;
    m_heading = data["initial"]["heading_angle[deg]"];
    m_heading = m_heading*pi/180.0;
    m_aileron = data["initial"]["aileron[deg]"];
    m_aileron = m_aileron*pi/180.0;
    m_elevator = data["initial"]["elevator[deg]"];
    m_elevator = m_elevator*pi/180.0;
    m_rudder = data["initial"]["rudder[deg]"];
    m_rudder = m_rudder*pi/180.0;
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
    
    // some calcs
    m_cw = m_wing_area/m_wing_span;
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
 
    // calculate V, alpha, beta, pbar, qbar, rbar
    double V = sqrt(pow(u,2) + pow(v,2) + pow(w,2));
    double alpha  = atan2(w, u);
    double beta = asin(v/V);
    double pbar = p*m_wing_span/(2*V);
    double qbar = q*m_cw/(2*V);
    double rbar = r*m_wing_span/(2*V);    

    // calculate CL, CD, and Cm
    double CL1 = m_CL0 + m_CL_a*alpha;
    double CL  = m_CL0 + m_CL_a*alpha + m_CL_qbar*qbar + m_CL_de*m_elevator;
    double CS  = m_CS_b*beta + m_CS_pbar*pbar + m_CS_rbar*rbar + m_CS_da*m_aileron + m_CS_dr*m_rudder;
    double CD  = m_CDL0 + m_CD_L*CL1 + m_CD_L2*CL1*CL1 + m_CD_S2*CS*CS + (m_CD_Lqbar*CL1 + m_CD_qbar)*qbar + (m_CD_Lde*CL1 + m_CD_de)*m_elevator + m_CD_de2*m_elevator*m_elevator;
    double Cl  = m_Cl_b*beta + m_Cl_pbar*pbar + (m_Cl_Lrbar*CL1 + m_Cl_rbar)*rbar + m_Cl_da*m_aileron + m_CL_de*m_elevator;
    double Cm  = m_Cm0 + m_Cm_a*alpha + m_Cm_qbar*qbar + m_Cm_de*m_elevator;
    double Cn  = m_Cn_b*beta + (m_Cn_Lpbar*CL1 + m_Cn_pbar)*pbar + m_Cn_rbar*rbar + (m_Cn_Lda*CL1 + m_Cn_da)*m_aileron + m_Cn_dr*m_rudder; 

    // get rho, mu, Re
    get_atmospheric_properties_english(-zf, m_atm);
    double rho = m_atm.density;
    double mu  = m_atm.dynamic_viscosity;
    //double Re  = rho*V*m_ref_length/mu;
    double ca = cos(alpha);
    double cb = cos(beta);
    double sa = sin(alpha);
    double sb = sin(beta);

    // update forces and moments
    ans[0] =  -0.5*rho*pow(V,2)*m_wing_area*(CD*ca*cb + CS*ca*sb - CL*sa); // F_xb
    ans[1] =   0.5*rho*pow(V,2)*m_wing_area*(CS*cb - CD*sb); // F_yb
    ans[2] =  -0.5*rho*pow(V,2)*m_wing_area*(CD*sa*cb + CS*sa*sb + CL*ca); // F_zb
    ans[3] =   0.5*rho*pow(V,2)*m_wing_area*m_wing_span*(Cl*ca*cb - Cm*ca*sb - Cn*sa); // 
    ans[4] =   0.5*rho*pow(V,2)*m_wing_area*m_wing_span*(Cl*sb + Cm*cb); // 
    ans[5] =   0.5*rho*pow(V,2)*m_wing_area*m_wing_span*(Cl*sa*cb - Cm*sa*sb + Cn*ca); // 
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
    
    
    int size = 13;
    double* y0 = new double[size];
    y0[0]  = sqrt((m_V*m_V - pow(m_V*sin(m_beta), 2))/(1 + pow(tan(m_alpha), 2)));         // u
    y0[1]  = m_V*sin(m_beta);              // v
    y0[2]  = y0[0]*tan(m_alpha);           // w
    y0[3]  = m_p;              // p
    y0[4]  = m_q;              // q
    y0[5]  = m_r;              // r
    y0[6]  = 0.0;              // xf
    y0[7]  = 0.0;              // yf
    y0[8]  = -m_altitude;      // zf
    y0[9]  = m_bank;           // phi
    y0[10] = m_elv_angle;      // theta
    y0[11] = m_heading;        // psi

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

    //print results in a file
    FILE* init_file = fopen("initial_state.txt", "w");
    fprintf(init_file, "Initial State\n");
    fprintf(init_file, "  Time[s]              u[ft/s]            v[ft/s]                w[ft/s]              p[rad/s]             q[rad/s]             r[rad/s]             x[ft]                y[ft]                z[ft]                e0                   ex                   ey                   ez\n");
    fprintf(init_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, y0[0], y0[1], y0[2],y0[3],y0[4], y0[5], y0[6], y0[7], y0[8], y0[9], y0[10], y0[11], y0[12]);
        
    // double* y = new double[13];
}

void aircraft::exercise_6_1()
{
    FILE* out_file = fopen("f_16.txt", "w");
    fprintf(out_file, "  Time[s]              u[ft/s]            v[ft/s]                w[ft/s]              p[rad/s]             q[rad/s]             r[rad/s]             x[ft]                y[ft]                z[ft]                e0                   ex                   ey                   ez\n");
    
    double t0 = 0.0;
    
    do {
        aircraft_rk4(t0, y0, m_time_step, size, y);
        
        fprintf(out_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, y0[0], y0[1], y0[2],y0[3],y0[4], y0[5], y0[6], y0[7], y0[8], y0[9], y0[10], y0[11], y0[12]);
        
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

