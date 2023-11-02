#include "aircraft.h"

// update aerodynamic coeifficint calcs in aerodynamics
// then update forces in func function


aircraft::aircraft(string filename)
{
    std::ifstream f(filename);
    m_input = json::parse(f);
    
    // read in json m_input and assign to existing class attributes
    printf("Reading Aircraft info from file:\n");
    
    // simulation
    m_time_step = m_input["simulation"]["time_step[s]"];
    m_total_time = m_input["simulation"]["total_time[s]"];
    
    // aircraft
    m_wing_area = m_input["aircraft"]["wing_area[ft^2]"];
    m_wing_span = m_input["aircraft"]["wing_span[ft]"];
    m_cw = m_wing_area/m_wing_span;
    m_weight = m_input["aircraft"]["weight[lbf]"];
    m_Ixx = m_input["aircraft"]["Ixx[slug-ft^2]"];
    m_Iyy = m_input["aircraft"]["Iyy[slug-ft^2]"];
    m_Izz = m_input["aircraft"]["Izz[slug-ft^2]"];
    m_Ixy = m_input["aircraft"]["Ixy[slug-ft^2]"];
    m_Ixz = m_input["aircraft"]["Ixz[slug-ft^2]"];
    m_Iyz = m_input["aircraft"]["Iyz[slug-ft^2]"];
    m_hx = m_input["aircraft"]["hx[slug-ft^2/s]"];
    m_hy = m_input["aircraft"]["hy[slug-ft^2/s]"];
    m_hz = m_input["aircraft"]["hz[slug-ft^2/s]"];

    json CG_shift = m_input["aircraft"]["CG_shift[ft]"];
    m_CG_shift[0] = CG_shift[0];
    m_CG_shift[1] = CG_shift[1];
    m_CG_shift[2] = CG_shift[2];    

    json thrust_location = m_input["aircraft"]["thrust"]["location[ft]"];
    m_thrust_loc[0] = thrust_location[0];
    m_thrust_loc[1] = thrust_location[1];
    m_thrust_loc[2] = thrust_location[2];

    json thrust_direction = m_input["aircraft"]["thrust"]["direction"];
    m_thrust_dir[0] = thrust_direction[0];
    m_thrust_dir[1] = thrust_direction[1];
    m_thrust_dir[2] = thrust_direction[2];
    vector_normalize_3(m_thrust_dir);

    m_T0 = m_input["aircraft"]["thrust"]["T0[lbf]"];
    m_T1 = m_input["aircraft"]["thrust"]["T1[lbf-s/ft]"];
    m_T2 = m_input["aircraft"]["thrust"]["T2[lbf-s^2/ft^2]"];
    m_a = m_input["aircraft"]["thrust"]["a"];

    

    // initial values for needed whether initialized with state or trim
    m_V = m_input["initial"]["airspeed[ft/s]"];
    m_altitude = m_input["initial"]["altitude[ft]"];
    m_heading = m_input["initial"]["heading[deg]"];
    m_heading *= pi/180.0;
 
    
    // some calcs
   
    Atmosphere atm;
    get_atmospheric_properties_english(0.0, atm);
    m_rho0 = atm.density; // [slug/ft^3]


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

     // aerodynamics
    m_CL0 = m_input["aerodynamics"]["CL"]["0"];
    m_CL_a = m_input["aerodynamics"]["CL"]["alpha"];
    m_CL_qbar = m_input["aerodynamics"]["CL"]["qbar"];
    m_CL_de = m_input["aerodynamics"]["CL"]["de"];

    m_CS_b = m_input["aerodynamics"]["CS"]["beta"];
    m_CS_pbar = m_input["aerodynamics"]["CS"]["pbar"];
    m_CS_Lpbar = m_input["aerodynamics"]["CS"]["Lpbar"];
    m_CS_rbar = m_input["aerodynamics"]["CS"]["rbar"];
    m_CS_da = m_input["aerodynamics"]["CS"]["da"];
    m_CS_dr = m_input["aerodynamics"]["CS"]["dr"];

    m_CDL0 = m_input["aerodynamics"]["CD"]["L0"];
    m_CD_L = m_input["aerodynamics"]["CD"]["L"];
    m_CD_L2 = m_input["aerodynamics"]["CD"]["L2"];
    m_CD_S2 = m_input["aerodynamics"]["CD"]["S2"];
    m_CD_qbar = m_input["aerodynamics"]["CD"]["qbar"];
    m_CD_Lqbar = m_input["aerodynamics"]["CD"]["Lqbar"];
    m_CD_de = m_input["aerodynamics"]["CD"]["de"];
    m_CD_Lde = m_input["aerodynamics"]["CD"]["Lde"];
    m_CD_de2 = m_input["aerodynamics"]["CD"]["de2"];
    
    m_Cl_b = m_input["aerodynamics"]["Cl"]["beta"];
    m_Cl_pbar = m_input["aerodynamics"]["Cl"]["pbar"];
    m_Cl_rbar = m_input["aerodynamics"]["Cl"]["rbar"];
    m_Cl_Lrbar = m_input["aerodynamics"]["Cl"]["Lrbar"];
    m_Cl_da = m_input["aerodynamics"]["Cl"]["da"];
    m_Cl_dr = m_input["aerodynamics"]["Cl"]["dr"];

    m_Cm0 = m_input["aerodynamics"]["Cm"]["0"];
    m_Cm_a = m_input["aerodynamics"]["Cm"]["alpha"];
    m_Cm_qbar = m_input["aerodynamics"]["Cm"]["qbar"];
    m_Cm_de = m_input["aerodynamics"]["Cm"]["de"];

    m_Cn_b = m_input["aerodynamics"]["Cn"]["beta"];
    m_Cn_pbar = m_input["aerodynamics"]["Cn"]["pbar"];
    m_Cn_Lpbar = m_input["aerodynamics"]["Cn"]["Lpbar"];
    m_Cn_rbar = m_input["aerodynamics"]["Cn"]["rbar"];
    m_Cn_da = m_input["aerodynamics"]["Cn"]["da"];
    m_Cn_Lda = m_input["aerodynamics"]["Cn"]["Lda"];
    m_Cn_dr = m_input["aerodynamics"]["Cn"]["dr"];

    m_size = 13;
    m_controls = new double[4];

    init_sim();  


    
}

void aircraft::init_sim()
{
    
    if (m_input["initial"]["type"] == "state"){
        init_from_state();
    }
    else {
        
        init_from_trim();
    }

}

void aircraft::init_from_state(){
    printf("Initializing from a given state");

    m_size = 13;
    m_elv_angle = m_input["initial"]["state"]["elevation_angle[deg]"];
    m_elv_angle = m_elv_angle*pi/180.0;
    m_bank_angle = m_input["initial"]["state"]["bank_angle[deg]"];
    m_bank_angle *= pi/180.0;
    m_alpha = m_input["initial"]["state"]["alpha[deg]"];
    m_alpha *= pi/180.0;
    m_beta = m_input["initial"]["state"]["beta[deg]"];
    m_beta *= pi/180.0;
    m_p = m_input["initial"]["state"]["p[deg/s]"];
    m_p *= pi/180.0;
    m_q = m_input["initial"]["state"]["q[deg/s]"];
    m_q *= pi/180.0;
    m_r = m_input["initial"]["state"]["r[deg/s]"];
    m_r *= pi/180.0;
    
    m_throttle = m_input["initial"]["state"]["throttle"];
    m_da = m_input["initial"]["state"]["aileron[deg]"];
    m_de = m_input["initial"]["state"]["elevator[deg]"];
    m_dr = m_input["initial"]["state"]["rudder[deg]"]; 


    m_controls = new double[4];
    m_controls[0] = m_da;
    m_controls[1] = m_de;
    m_controls[2] = m_dr;
    m_controls[3] = m_throttle;    
    
    m_initial_state = new double[m_size];      
    m_initial_state[0]  = m_V*cos(m_alpha)*cos(m_beta);  // u
    m_initial_state[1]  = m_V*sin(m_beta);               // v
    m_initial_state[2]  = m_V*sin(m_alpha)*cos(m_beta);  // w
    m_initial_state[3]  = m_p;              // p
    m_initial_state[4]  = m_q;              // q
    m_initial_state[5]  = m_r;              // r
    m_initial_state[6]  = 0.0;              // xf
    m_initial_state[7]  = 0.0;              // yf
    m_initial_state[8]  = -m_altitude;      // zf
    m_initial_state[9]  = m_bank_angle;           // phi
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

void aircraft::init_from_trim(){
    double sa,ca, sb, cb, sp, cp, sg, st, ct, u, v, w, pqr_constant, var, grav;
    double term, bottom, big_term, theta1, theta2, check1, check2;
    double Rmax, n_a, qbar, CLift, thrust_mag;
    double G[6], delta_G[6] ;
    double R[6], R_up[6], R_down[6], R_copy[6], FM[6];

    
    
    // read in json values
    m_trim_type = m_input["initial"]["trim"]["type"];
    m_finite_diff_step = m_input["initial"]["trim"]["solver"]["finite_difference_step_size"];
    m_relaxation = m_input["initial"]["trim"]["solver"]["relaxation_factor"];
    m_tol = m_input["initial"]["trim"]["solver"]["tolerance"];
    m_max_iter = m_input["initial"]["trim"]["solver"]["max_iterations"];
    m_verbose = m_input["initial"]["trim"]["solver"]["verbose"];
    cout<< "Verbose = " << m_verbose << endl;
    

    
    get_atmospheric_properties_english(m_altitude, m_atm);
    double Rho = m_atm.density;
    //double CW = 0.5*Rho*m_V*m_V*m_wing_area*m_weight;

    bool converged = false;
    
    // step 1: initialize alpha, beta, da, de, dr, tau = 0
    double alpha = 0.0;    
    double beta = 0.0;
    
    double da = 0.0;    // da
    double de = 0.0;    // de
    double dr = 0.0;    // dr
    double tau = 0.0;    // tau

    

    double* y = new double[9];
    y[0] = m_V*cos(alpha)*cos(beta);
    y[1] = m_V*sin(beta);
    y[2] = m_V*sin(alpha)*cos(beta);

    y[3] = 0.0;
    y[4] = 0.0;
    y[5] = 0.0;

    y[6] = 0.0; 
    y[7] = 0.0; 
    y[8] = -m_altitude; 


    //initialize angles at zero
    double gamma = 0.0;
    double theta = 0.0;
    double phi = 0.0;

    bool solve_Elevation = false;
    // check to see if given elevation angle
    if (key_exists(m_input["initial"]["trim"],"elevation_angle[deg]")){
        theta = m_input["initial"]["trim"]["elevation_angle[deg]"];
        theta = theta*pi/180.0;
    }else{
    // if given climb angle, conversion to elevation will be needed
        gamma = m_input["initial"]["trim"]["climb_angle[deg]"];
        gamma = gamma*pi/180.0;
        solve_Elevation = true;
    }

    bool solve_traditional = true;
    // check to see if given phi
    if (key_exists(m_input["initial"]["trim"],"bank_angle[deg]")){
        phi = m_input["initial"]["trim"]["bank_angle[deg]"];
        phi *= pi/180.0;
    }else{
    // if given sideslip, conversion to phi will be needed
        beta = m_input["initial"]["trim"]["sideslip_angle[deg]"];
        beta = beta*pi/180.0;
        solve_traditional = false;
    }
    


    if (m_trim_type == "sct"){
        printf("\nTriming aircraft for a Steady Coorinated Turn\n");
        // requires bank angle and requires elv or climb angle
        if (m_verbose == true){
            printf("\nInitial phi   = %f\n", phi*180.0/pi);
            printf("Initial theta = %f\n\n", theta*180.0/pi);
        }

    }else if (m_trim_type == "shss"){
        printf("\nTrimming Aircraft for a Steady Heading Side Slip\n");
        // requires bank or side slip angle
        // requires elv or climb angle


        if (m_verbose == true){
            printf("\nInitial phi   = %f\n", phi*180.0/pi);
            printf("Initial theta = %f\n\n", theta*180.0/pi);
            printf("p[deg/s]        = %20.12e\n", y[3]*180.0/pi);
            printf("q[deg/s]        = %20.12e\n", y[4]*180.0/pi);
            printf("r[deg/s]        = %20.12e\n\n", y[5]*180.0/pi);
        }
        
    }


    // build G vector
    G[0] = alpha;

    if (solve_traditional == true){
        G[1] = beta;
        var = phi;
    }else{
        G[1] = phi;
        var = beta;
    }

    G[2] = da; 
    G[3] = de;    
    G[4] = dr; 
    G[5] = tau; 

        
    int iter = 0;
    do {
        //gamma = asin((u*st - (v*sp + w*cp)*ct)/m_V); 
        sg = sin(gamma);

        if (solve_traditional == true){
            beta = G[1];
            phi = var;
        }else{
            phi = G[1];
            beta = var;
        }
        
        alpha = G[0];
        
        sa = sin(alpha);
        ca = cos(alpha);
        sb = sin(beta);
        cb = cos(beta);
        sp = sin(phi);
        cp = cos(phi);


        y[0] = m_V*ca*cb; // u
        y[1] = m_V*sb;    // v
        y[2] = m_V*sa*cb; // w

        

        u = y[0]; // u
        v = y[1];    // v
        w = y[2]; // w
        
        // update theta
        if (solve_Elevation == true){
            printf("\nSolving for theta given a climb angle\n");
            term = (v*sp + w*cp);
            bottom = u*u + term*term;
            big_term = term*sqrt(u*u + term*term - m_V*m_V*sg*sg);
            theta1 = asin((u*m_V*sg + big_term)/bottom);
            theta2 = asin((u*m_V*sg - big_term)/bottom);
            check1 = m_V*sg -u*sin(theta1) + term*cos(theta1);
            check2 = m_V*sg -u*sin(theta2) + term*cos(theta2);
            if (abs(check1) < abs(check2)){ theta = theta1;}
            else{theta = theta2;}

            if (m_verbose == true){
                printf("theta 1        = %20.12e\n", theta1*180.0/pi);
                printf("theta 2        = %20.12e\n", theta2*180.0/pi);
                printf("correct theta  = %20.12e\n\n", theta*180.0/pi);}   
        }
        


        if (m_trim_type == "sct"){
        
            grav = gravity_english(m_altitude);// 
            st = sin(theta);
            ct = cos(theta);
            pqr_constant = grav*sp*ct/(y[0]*ct*cp + y[2]*st);
            y[3] = -pqr_constant*st;    // p
            y[4] = pqr_constant*sp*ct;  // q
            y[5] = pqr_constant*cp*ct;  // r

            if (m_verbose == true){
                printf("\nUpdating p, q and r\n");
                printf("p[deg/s]        = %20.12e\n", y[3]*180.0/pi);
                printf("q[deg/s]        = %20.12e\n", y[4]*180.0/pi);
                printf("r[deg/s]        = %20.12e\n\n", y[5]*180.0/pi);}
        }

        // get R unperturbed
        if (m_verbose == true){
            if (solve_traditional == true){
                ("\n\nG  = [alpha, beta, da, de, dr, tau]\n");}
            else{
                ("\n\nG  = [alpha, phi, da, de, dr, tau]\n");
            }
        }
        
        calc_R(G, y, var, theta, R, solve_traditional);
        
        // initialize jacobian
        if (m_verbose == true){
            printf("Building Jacobian Matrix:\n");}
        
        double** J = new double* [6];
        for (int i = 0; i< 6;i++){
            J[i] = new double[6];}
        for (int i = 0; i< 6;i++){
            for (int j = 0; j< 6;j++){
                J[i][j] = 0.0;}}

        if (m_verbose == true){
            printf("Finite Difference step size = %f \n\n", m_finite_diff_step);}
        
        // loop through perturbing G elements
        for (int i = 0; i < 6; i++){
            if (m_verbose == true){
                printf("-------------------------------------\n");
                printf("Computing Gradient relative to G[%d]\n",i);
                printf("-------------------------------------\n\n");}
            
            // step up
            G[i] += m_finite_diff_step;

            // calc R_up
            if (m_verbose == true){
                printf("Positive Finite Difference Step\n");}
            
            calc_R(G, y, var, theta, R_up, solve_traditional);
            
            // step down
            G[i] -= 2.0*m_finite_diff_step;

            // calc R_down
            if (m_verbose == true){
                printf("Negative Finite Difference Step\n");}
            
            calc_R(G, y, var, theta, R_down, solve_traditional);

            // update column
            for (int j = 0; j < 6; j++){
                J[j][i] = (R_up[j] - R_down[j])/(2*m_finite_diff_step);
            }

            // put G back to normal
            G[i] += m_finite_diff_step;
        }
    
        if (m_verbose == true){
                printf("Jacobian Matrix:\n");
                for (int i = 0; i< 6;i++){
                    array_print(J[i], 6);}}

                    

        array_copy(R,R_copy,6);
        // solve for delta G
        matrix_AxB_solve(J, R_copy, 6, delta_G);

        for (int i = 0;i<6;i++){
            delta_G[i]= -delta_G[i];}

        if (m_verbose == true){
            printf("\nDelta G = ");array_print(delta_G, 6);}

        // step G with relaxation factor
        for (int i = 0; i< 6;i++){
            G[i] += delta_G[i]*m_relaxation;}
        
        if (m_verbose == true){
            printf("Relaxation Factor = %f\n", m_relaxation);
            printf("New G  = ");array_print(G, 6);printf("\n");}

        // calc CL 
        //qbar = y[4]*m_cw/(2.0*m_V); 
        //CLift  = m_CL0 + m_CL_a*G[0] + m_CL_qbar*qbar + m_CL_de*G[3];
        
        // calc load factor eq 16.32
        aerodynamics_aircraft(y, FM);
        n_a = (-FM[2]*cos(G[0]) + FM[0]*sin(G[0]))/m_weight;

        // calc thrust
        thrust_mag = G[5]*pow(Rho/m_rho0, m_a)*(m_T0 + m_T1*m_V + m_T2*m_V*m_V);
        
        Rmax = 0.0;
        for (int i = 0; i< 6;i++){
                Rmax = max(Rmax, abs(R[i]));}

        if (Rmax < m_tol){
            converged = true;}

        iter += 1;

        if (iter > m_max_iter){
            converged = true;}
        

        if (m_verbose == true){
            printf("Iteration     Throttle              Alpha[deg]         Beta[deg]  \n");
            cout<< iter << "             " << G[5]<< "      " <<  G[0]*180.0/pi << "      " << G[1]*180.0/pi  << endl;
            
            printf("\nAileron[deg]          Elevator[deg]           Rudder[deg]\n");
            cout<< G[2]*180.0/pi<< "      " << G[3]*180.0/pi <<  "      " << G[4]*180.0/pi << endl;
            
            printf("\np[deg/s]               q[deg/s]             r[deg/s]                 phi[deg]                theta[deg]        \n");      
            cout<< y[3]*180.0/pi << "      " << y[4]*180.0/pi << "      " << y[5]*180.0/pi << "      " << phi*180.0/pi << "      " << theta*180.0/pi << endl;
            
            printf("\nLoad Factor           Max Residual\n");
            cout<< n_a << "             " <<   Rmax  <<  endl;}

        delete[] J;

    } while (converged == false);//(iter<2);//

    // calc thrust
    // write to init state vector

    
    printf("\n\n--------------- Trim Solution ---------------\n");
    printf("elevation angle[deg] = %20.12e\n", theta*180.0/pi);
    printf("bank angle[deg]      = %20.12e\n", phi*180.0/pi);
    printf("alpha[deg]           = %20.12e\n", G[0]*180.0/pi);
    printf("beta[deg]            = %20.12e\n", beta*180.0/pi);
    printf("p[deg/s]             = %20.12e\n", y[3]*180.0/pi);
    printf("q[deg/s]             = %20.12e\n", y[4]*180.0/pi);
    printf("r[deg/s]             = %20.12e\n", y[5]*180.0/pi);
    printf("aileron[deg]         = %20.12e\n", G[2]*180.0/pi);
    printf("elevator[deg]        = %20.12e\n", G[3]*180.0/pi);
    printf("rudder[deg]          = %20.12e\n", G[4]*180.0/pi);
    printf("Trottle              = %20.12e\n", G[5]);
    printf("Thrust[lbf]          = %20.12e\n\n", thrust_mag);
    printf("check");

        

}

void aircraft::calc_R(double G[6], double* y, double var, double theta, double ans[6], bool solve_traditional){
    
    double alpha, beta, phi, da, de, dr, tau, sp, cp, st, ct, grav, pqr_constant;
    double* FM = new double[6];
    
    if (solve_traditional == true){
        // if traditional, phi is passsed in as var, beta is G[1]
        phi = var;
        
        // update alpha, beta, controls
        beta  = G[1];
      
        
    }else{
        // if non traditional, beta is passed in as var and phi is passed in as G[1]
        beta = var;
        
        phi   = G[1];

    }
        
        alpha = G[0];
        m_controls[0]    = G[2]; // da
        m_controls[1]    = G[3]; // de
        m_controls[2]    = G[4]; // dr
        m_controls[3]    = G[5]; // dr

        sp = sin(phi);
        cp = cos(phi);
        st = sin(theta);
        ct = cos(theta);
        grav = gravity_english(m_altitude);
    
    // step 4: Calculate the body-fixed velocities from Eq. (14.9) for the traditional definition of sideslip
    y[0] = m_V*cos(alpha)*cos(beta); // u
    y[1] = m_V*sin(beta);            // v
    y[2] = m_V*sin(alpha)*cos(beta); // w
    
    // pqr_constant = grav*sp*ct/(y[0]*ct*cp + y[2]*st);
    // y[3] = -pqr_constant*st;    // p
    // y[4] = pqr_constant*sp*ct;  // q
    // y[5] = pqr_constant*cp*ct;  // r

    // call calc aero
    aerodynamics_aircraft(y, FM);

    double u = y[0];
    double v = y[1];
    double w = y[2];
    double p = y[3];
    double q = y[4];
    double r = y[5];
    

    //put FM in the following
    ans[0] = FM[0] - m_weight*st    + (r*v - q*w)*m_weight/grav;
    ans[1] = FM[1] + m_weight*sp*ct + (p*w - r*u)*m_weight/grav;
    ans[2] = FM[2] + m_weight*cp*ct + (q*u - p*v)*m_weight/grav;;
    ans[3] = FM[3] - m_hz*q + m_hy*r + (m_Iyy - m_Izz)*q*r + m_Iyz*(q*q - r*r) + m_Ixz*p*q - m_Ixy*p*r;
    ans[4] = FM[4] + m_hz*p - m_hx*r + (m_Izz - m_Ixx)*p*r + m_Ixz*(r*r - p*p) + m_Ixy*q*r - m_Iyz*p*q;
    ans[5] = FM[5] - m_hy*p + m_hx*q + (m_Ixx - m_Iyy)*p*q + m_Ixy*(p*p - q*q) + m_Iyz*p*r - m_Ixz*q*r;
    
    if (m_verbose == true){
        printf("G  = ");array_print(G, 6);;
        printf("y  = ");array_print(y, 9);;
        printf("FM = ");array_print(FM, 6);;
        printf("R  = ");array_print(ans, 6);printf("\n");}



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

    
    double da = m_controls[0];
    double de = m_controls[1];
    double dr = m_controls[2];   
    double tau = m_controls[3]; // throttle setting

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

    double constant = 0.5*rho*V*V*m_wing_area;

    // update forces and moments
    ans[0] =  -constant*(CD*ca*cb + CS*ca*sb - CL*sa); // F_xb
    ans[1] =   constant*(CS*cb - CD*sb); // F_yb
    ans[2] =  -constant*(CD*sa*cb + CS*sa*sb + CL*ca); // F_zb
    ans[3] =   constant*m_wing_span*Cl; // 
    ans[4] =   constant*m_cw*Cm; // 
    ans[5] =   constant*m_wing_span*Cn; // 

    // update and add thrust with thrust vector
    
    m_thrust_mag = tau*pow(rho/m_rho0, m_a)*(m_T0 + m_T1*V + m_T2*V*V);
    double thrust[3];

    thrust[0] = m_thrust_dir[0]*m_thrust_mag;
    thrust[1] = m_thrust_dir[1]*m_thrust_mag;
    thrust[2] = m_thrust_dir[2]*m_thrust_mag;

    ans[0] += thrust[0];
    ans[1] += thrust[1];
    ans[2] += thrust[2];

    // get propulsive moment
    double M_p[3];
    vector_cross_3(m_thrust_loc, thrust, M_p);

    ans[3] += M_p[0];
    ans[4] += M_p[1];
    ans[5] += M_p[2]; 
    
    // update moments with CG shift
    double shift[3];
    vector_cross_3(m_CG_shift, &ans[0], shift);
    
    ans[3] -= shift[0];
    ans[4] -= shift[1];
    ans[5] -= shift[2];

}

void aircraft::aircraft_rk4_func(double t, double* y, double* ans)
{

    double* FM = new double[6];
    // update aerodynamic m_input
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
    aircraft_rk4_func(t0, y0, m_k1);
    fprintf(m_check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, m_k1[0], m_k1[1], m_k1[2],m_k1[3],m_k1[4], m_k1[5], m_k1[6], m_k1[7], m_k1[8], m_k1[9], m_k1[10], m_k1[11], m_k1[12]);

    for (int i=0; i<size; i++) {ans[i] = y0[i] + 0.5*dt*m_k1[i];}
    
    aircraft_rk4_func(t0+0.5*dt, ans, m_k2);
    fprintf(m_check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, m_k2[0], m_k2[1], m_k2[2],m_k2[3],m_k2[4], m_k2[5], m_k2[6], m_k2[7], m_k2[8], m_k2[9], m_k2[10], m_k2[11], m_k2[12]);

    for (int i=0; i<size; i++) {ans[i] = y0[i] + 0.5*dt*m_k2[i];}

    aircraft_rk4_func(t0+0.5*dt, ans, m_k3);
    fprintf(m_check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, m_k3[0], m_k3[1], m_k3[2],m_k3[3],m_k3[4], m_k3[5], m_k3[6], m_k3[7], m_k3[8], m_k3[9], m_k3[10], m_k3[11], m_k3[12]);

    for (int i=0; i<size; i++) {ans[i] = y0[i] + dt*m_k3[i];}

    aircraft_rk4_func(t0+dt, ans, m_k4);
    fprintf(m_check_file, "%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", t0, m_k4[0], m_k4[1], m_k4[2],m_k4[3],m_k4[4], m_k4[5], m_k4[6], m_k4[7], m_k4[8], m_k4[9], m_k4[10], m_k4[11], m_k4[12]);

    for (int i=0; i<size; i++){ans[i] = y0[i] + dt/6.0*(m_k1[i] + 2.0*m_k2[i] + 2.0*m_k3[i] + m_k4[i]);}
    
    fprintf(m_check_file, "\n"); 
}

void aircraft::run_sim()
{
    m_check_file = fopen("check_rk4.txt", "w");
    fprintf(m_check_file, "  Time[s]              udot                vdot                 wdot                  pdot               qdot                 rdot                 xdot                  ydot                 zdot                   e0dot                exdot               eydot              ezdot\n");
    
    FILE* out_file = fopen("hoch_7_6.txt", "w");
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

