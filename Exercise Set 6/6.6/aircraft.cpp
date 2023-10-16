#include "aircraft.h"

// update aerodynamic coeifficint calcs in aerodynamics
// then update forces in func function


aircraft::aircraft(string filename)
{
    std::ifstream f(filename);
    json data = json::parse(f);
    
    // read in json data and assign to existing class attributes
    
    // simulation
    m_time_step = data["simulation"]["time_step[s]"];
    m_total_time = data["simulation"]["total_time[s]"];

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
    json CG_shift = data["aircraft"]["CG_shift[ft]"];
    m_CG_shiftx = CG_shift[0];
    m_CG_shifty = CG_shift[1];
    m_CG_shiftz = CG_shift[2];
    json thrust_location = data["aircraft"]["thrust"]["location[ft]"];
    m_thrust_locx = thrust_location[0];
    m_thrust_locy = thrust_location[1];
    m_thrust_locz = thrust_location[2];
    json thrust_direction = data["aircraft"]["thrust"]["direction"];
    m_thrust_dirx = thrust_direction[0];
    m_thrust_diry = thrust_direction[1];
    m_thrust_dirz = thrust_direction[2];
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
    m_bank *= pi/180.0;
    m_alpha = data["initial"]["alpha[deg]"];
    m_alpha *= pi/180.0;
    m_beta = data["initial"]["beta[deg]"];
    m_beta *= pi/180.0;
    m_p = data["initial"]["p[deg/s]"];
    m_p *= pi/180.0;
    m_q = data["initial"]["q[deg/s]"];
    m_q *= pi/180.0;
    m_r = data["initial"]["r[deg/s]"];
    m_r *= pi/180.0;
    m_heading = data["initial"]["heading_angle[deg]"];
    m_heading *= pi/180.0;
    m_da = data["initial"]["aileron[deg]"];
    m_da *= pi/180.0;
    m_de = data["initial"]["elevator[deg]"];
    m_de *= pi/180.0;
    m_dr = data["initial"]["rudder[deg]"];
    m_dr *= pi/180.0;
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
    Atmosphere atm;
    get_atmospheric_properties_english(0.0, atm);
    m_rho0 = atm.density; // [slug/ft^3]
    //cout << "rho0" << endl;
    //cout << m_rho0 << endl;

    // make I matrix and invert it
    m_I[0][0] = m_Ixx;
    m_I[0][1] = -m_Ixy;
    m_I[0][2] = -m_Ixz;
    m_I[1][0] = -m_Ixy;
    m_I[1][1] = m_Iyy;
    m_I[1][2] = -m_Iyz;
    m_I[2][0] = -m_Ixz;
    m_I[2][1] = -m_Iyz;
    m_I[2][2] = m_Izz;

    matrix_invert_3x3(m_I, m_I_inv);
    //cout << "I" << endl;
    //array_print_3x3(m_I);
    
    //cout << "I inv" << endl;
    //array_print_3x3(m_I_inv);


    m_h[0][0] = 0.0;
    m_h[0][1] = -m_hz;
    m_h[0][2] = m_hy;
    m_h[1][0] = m_hz;
    m_h[1][1] = 0.0;
    m_h[1][2] = -m_hx;
    m_h[2][0] = -m_hy;
    m_h[2][1] = m_hx;
    m_h[2][2] = 0.0;

    //cout << "h" << endl;
    //array_print_3x3(m_h);
    
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
 
    // calculate V, alpha, beta, pbar, qbar, rbar
    double V = sqrt(u*u + v*v + w*w);
    double alpha  = atan2(w, u);
    double beta = asin(v/V);
    double pbar = p*m_wing_span/(2.0*V);
    double qbar = q*m_cw/(2.0*V);
    double rbar = r*m_wing_span/(2.0*V); 

    double tau = m_controls[0]; // throttle setting
    double da = m_controls[1];
    double de = m_controls[2];
    double dr = m_controls[3];   

    // calculate CL, CD, and Cm
    double CL1 = m_CL0 + m_CL_a*alpha;
    double CL  = CL1 + m_CL_qbar*qbar + m_CL_de*de;
    double CS  = m_CS_b*beta + m_CS_pbar*pbar + m_CS_rbar*rbar + m_CS_da*da + m_CS_dr*dr;
    double CD  = m_CDL0 + m_CD_L*CL1 + m_CD_L2*CL1*CL1 + m_CD_S2*CS*CS + (m_CD_Lqbar*CL1 + m_CD_qbar)*qbar + (m_CD_Lde*CL1 + m_CD_de)*de + m_CD_de2*de*de;
    double Cl  = m_Cl_b*beta + m_Cl_pbar*pbar + (m_Cl_Lrbar*CL1 + m_Cl_rbar)*rbar + m_Cl_da*da + m_Cl_dr*dr;
    double Cm  = m_Cm0 + m_Cm_a*alpha + m_Cm_qbar*qbar + m_Cm_de*de;
    double Cn  = m_Cn_b*beta + (m_Cn_Lpbar*CL1 + m_Cn_pbar)*pbar + m_Cn_rbar*rbar + (m_Cn_Lda*CL1 + m_Cn_da)*da + m_Cn_dr*dr; 

    //cout << "de" << endl;
    //cout << de << endl;
   


    // get rho, mu, Re
    get_atmospheric_properties_english(-zf, m_atm);
    double rho = m_atm.density;
    //double Re  = rho*V*m_ref_length/mu;
    double ca = cos(alpha);
    double cb = cos(beta);
    double sa = sin(alpha);
    double sb = sin(beta);

    // update forces and moments
    ans[0] =  -0.5*rho*V*V*m_wing_area*(CD*ca*cb + CS*ca*sb - CL*sa); // F_xb
    ans[1] =   0.5*rho*V*V*m_wing_area*(CS*cb - CD*sb); // F_yb
    ans[2] =  -0.5*rho*V*V*m_wing_area*(CD*sa*cb + CS*sa*sb + CL*ca); // F_zb
    ans[3] =   0.5*rho*V*V*m_wing_area*m_wing_span*Cl; // 
    ans[4] =   0.5*rho*V*V*m_wing_area*m_cw*Cm; // 
    ans[5] =   0.5*rho*V*V*m_wing_area*m_wing_span*Cn; // 

    // update and add thrust with thrust vector
    double thrust[3];
    thrust[0] = tau*pow(rho/m_rho0, m_a)*(m_T0 + m_T1*V + m_T2*V*V);
    ans[0] += thrust[0];
    ans[1] += thrust[1];
    ans[2] += thrust[2];
    
    // update moments with CG shift
    double shift[3];
    vector_cross_3(m_CG_shift, &ans[0], shift);
    ans[3] -= shift[0];
    ans[4] -= shift[1];
    ans[5] -= shift[2];

    //cout<< "ans=" << endl;
    //array_print(ans,13);


    //calculuate proptulsive thrust
    // calculate location
    // cross thrust location by force p

    // get cg shift = moments - (new cg loc x forcexyz)
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

    // create pqr matrix
    double pqr[3];//] = new double[3];
    pqr[0] = p;
    pqr[1] = q;
    pqr[2] = r;
    //cout << "pqr" << endl;
    //array_print(pqr, 3);

    // multiply h matrix by pqr
    double* hpqr =  new double[3];
    matrix_vector_mult_3(m_h, pqr, hpqr);

    //cout << "hpqr" << endl;
    //array_print(hpqr, 3);

    double pqr_dot_stuff[3];
    pqr_dot_stuff[0] = hpqr[0] + Mxb + (m_Iyy - m_Izz)*q*r + m_Iyz*(q*q - r*r) + m_Ixz*p*q - m_Ixy*p*r;
    pqr_dot_stuff[1] = hpqr[1] + Myb + (m_Izz - m_Ixx)*p*r + m_Ixz*(r*r - p*p) + m_Ixy*q*r - m_Iyz*p*q;
    pqr_dot_stuff[2] = hpqr[2] + Mzb + (m_Ixx - m_Iyy)*p*q + m_Ixy*(p*p - q*q) + m_Iyz*p*r - m_Ixz*q*r;

    // multiply I inv by hpqr
    double* pqr_dot = new double[3];
    matrix_vector_mult_3(m_I_inv, pqr_dot_stuff, pqr_dot);
    
    // cout << "pqr_dot_stuff[0]" << endl;
    // cout << pqr_dot_stuff[0] << endl;
    // cout << "pqr_dot_stuff[1]" << endl;
    // cout << pqr_dot_stuff[1] << endl;
    // cout << "pqr_dot_stuff[2]" << endl;
    // cout << pqr_dot_stuff[2] << endl;
    

   

    ans[0]  = (g*Fxb/m_weight) + (g*2.0*(ex*ez - ey*e0)) + (r*v) - (q*w)  ; // udot
    ans[1]  = (g*Fyb/m_weight) + (g*2.0*(ey*ez + ex*e0)) + (p*w) - (r*u); // vdot
    ans[2]  = (g*Fzb/m_weight) + (g*(ez*ez + e0*e0 - ex*ex - ey*ey)) + (q*u) - (p*v); // wdot
    ans[3]  = pqr_dot[0]; // pdot
    ans[4]  = pqr_dot[1]; // qdot
    ans[5]  = pqr_dot[2]; //  rdot
    //cout << "pqr_dot" << endl;
    //array_print(pqr_dot, 3);

     
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
    
    double* dy = new double[size];
    // k1:
    aircraft_rk4_func(t0, y0, dy);
    fprintf(m_check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11], dy[12]);
    //array_print(dy, 13);
    // multiply k1 by dt 
    for (int i = 0; i < size; i++)
    {
        m_k1[i] = dt * dy[i];
        // use k1 to update y_temporary
        m_y_temp[i] = y0[i] + 0.5 * m_k1[i];
    }

    // k2:  
    aircraft_rk4_func(t0 + (0.5 * dt), m_y_temp, dy);
    fprintf(m_check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11], dy[12]);
    //array_print(dy, 13);
    // multiply k2 by dt and then update y_temp
    for (int i = 0; i < size; i++)
    {
        m_k2[i] = dt * dy[i];
        m_y_temp[i] = y0[i] + (0.5 * m_k2[i]);
    }
    
    // k3:
    aircraft_rk4_func(t0 + (0.5 * dt), m_y_temp, dy);
    fprintf(m_check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11], dy[12]);//fprintf(check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11]);
    //array_print(dy, 13);
    // multiply k2 by dt and then update y_temp
    for (int i = 0; i < size; i++)
    {  
        m_k3[i] = dt * dy[i];
        m_y_temp[i] = y0[i] + m_k3[i];
    }
    
    // k4:
    aircraft_rk4_func(t0 + dt, m_y_temp, dy);
    fprintf(m_check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11], dy[12]);//fprintf(check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, dy[0], dy[1], dy[2],dy[3],dy[4], dy[5], dy[6], dy[7], dy[8], dy[9], dy[10], dy[11]);
    fprintf(m_check_file, "\n"); 
    //array_print(dy, 13);
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
    
    m_controls = new double[4];
    m_controls[0] = m_throttle;
    m_controls[1] = m_da;
    m_controls[2] = m_de;
    m_controls[3] = m_dr;


    
    m_CG_shift[0] = m_CG_shiftx;
    m_CG_shift[1] = m_CG_shifty;
    m_CG_shift[2] = m_CG_shiftz;
    
    m_thrust_loc = new double[3];
    m_thrust_loc[0] = m_thrust_locx;
    m_thrust_loc[1] = m_thrust_locy;
    m_thrust_loc[2] = m_thrust_locz;
    
    m_thrust_dir = new double[3];
    m_thrust_dir[0] = m_thrust_dirx;
    m_thrust_dir[1] = m_thrust_diry;
    m_thrust_dir[2] = m_thrust_dirz;

    vector_normalize_3(m_thrust_dir);

    
    m_size = 13;
    m_initial_state = new double[m_size];
    m_initial_state[0]  = sqrt((m_V*m_V - pow(m_V*sin(m_beta), 2))/(1 + pow(tan(m_alpha), 2)));         // u
    //cout << " check init " << endl;
    //cout << setprecision(12) << m_initial_state[0] << endl;
    m_initial_state[1]  = m_V*sin(m_beta);              // v
    m_initial_state[2]  = m_initial_state[0]*tan(m_alpha);           // w
    m_initial_state[3]  = m_p;              // p
    m_initial_state[4]  = m_q;              // q
    m_initial_state[5]  = m_r;              // r
    m_initial_state[6]  = 0.0;              // xf
    m_initial_state[7]  = 0.0;              // yf
    m_initial_state[8]  = -m_altitude;      // zf
    m_initial_state[9]  = m_bank;           // phi
    m_initial_state[10] = m_elv_angle;      // theta
    m_initial_state[11] = m_heading;        // psi

    // convert phi, theta, and psi to a quat
    double* quat = new double[4];
    euler_to_quat(&m_initial_state[9], quat);

    // normalize the quat
    quat_norm(quat);

    // store the quat in y0
    m_initial_state[9] = quat[0];
    m_initial_state[10] = quat[1];
    m_initial_state[11] = quat[2];
    m_initial_state[12] = quat[3];

    double t0 = 0.0;
    //print results in a file
    FILE* init_file = fopen("initial_state.txt", "w");
    fprintf(init_file, "Initial State\n");
    fprintf(init_file, "  Time[s]              u[ft/s]            v[ft/s]                w[ft/s]              p[rad/s]             q[rad/s]             r[rad/s]             x[ft]                y[ft]                z[ft]                e0                   ex                   ey                   ez\n");
    fprintf(init_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, m_initial_state[0], m_initial_state[1], m_initial_state[2],m_initial_state[3],m_initial_state[4], m_initial_state[5], m_initial_state[6], m_initial_state[7], m_initial_state[8], m_initial_state[9], m_initial_state[10], m_initial_state[11], m_initial_state[12]);
        
}

void aircraft::exercise_6_6()
{
    m_check_file = fopen("check_rk4.txt", "w");
    fprintf(m_check_file, "  Time[s]              udot                vdot                 wdot                  pdot               qdot                 rdot                 xdot                  ydot                 zdot                   e0dot                exdot               eydot              ezdot\n");
    
    FILE* out_file = fopen("hoch_6_6.txt", "w");
    fprintf(out_file, "  Time[s]              u[ft/s]            v[ft/s]                w[ft/s]              p[rad/s]             q[rad/s]             r[rad/s]             x[ft]                y[ft]                z[ft]                e0                   ex                   ey                   ez\n");
    
    double t0 = 0.0;
    //cout << "check in exercise" << endl;
    //cout << setprecision(12) << m_initial_state[0] << endl;

    double* y0 = new double[m_size];
    array_copy(m_initial_state, y0, m_size);

    double* y = new double[m_size];
    do {
        aircraft_rk4(t0, y0, m_time_step, m_size, y);
        
        
        fprintf(out_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, y0[0], y0[1], y0[2],y0[3],y0[4], y0[5], y0[6], y0[7], y0[8], y0[9], y0[10], y0[11], y0[12]);
        
        // y = y0 ( copy y0 into y)
        array_copy(y, y0, m_size);

        // re-normalize the quat
        quat_norm(&y[9]);
        
        // add time step
        t0 += m_time_step;

    } while (t0 <= m_total_time); //0.1); //
    //while (t0<0.005);
    
    fclose(m_check_file);
    fclose(out_file);
    
    
}

