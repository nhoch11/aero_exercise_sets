#ifndef CANNONBALL_H
#define CANNONBALL_H


#include "hoch_helper.h"
#include <string>


class cannonball
{
public: // private variables and function can only be accessed inside the class
    cannonball(std::string filename); // class constructor
    ~cannonball(){}; // class destructor (default: releases memory)
    
    void exercise_4_3();
    
private:
    
    // declare a struct from hoch_helper
    Atmosphere m_atm;

    // declare json variables
    double m_time_step;
    double m_ref_area;
    double m_ref_length;
    double m_init_V;
    double m_init_altitude;
    double m_init_theta;
    double m_p_initial;
    double m_q_initial;
    double m_r_initial;
    double m_weight;
    double m_Ixx;
    double m_Iyy;
    // aerodynamic coefficients
    double m_CLa;
    double m_CD0;
    double m_CD2;
    double m_Cma;
    double m_Cmq;
    double m_Clp;
    double m_Cl0;

    double m_Re_points[44] = {38.3119, 46.8488, 61.8966, 83.0529, 113.1784, 159.0790, 253.0613, 
    415.2212, 736.0977, 1152.9964, 1951.2934, 3623.5961, 5945.5707, 9605.6527, 15760.8866, 
    28376.3799, 42431.4792, 56934.6092, 70706.9926, 89180.2770, 105728.8597, 114234.0828, 
    117824.7587, 123423.4976, 137543.7498, 165609.7957, 193326.0255, 222215.4207, 267558.8711, 
    312337.1592, 359010.8111, 419094.3712, 496862.8878, 589062.3835, 742963.9508, 951688.0666, 
    1276976.1784, 1939252.7171, 2725736.5074, 3831186.8496, 6094614.4435, 18570226.6481, 
    27768242.0987, 10475175.0779};

    double m_cd_points[44] = {1.751558828, 1.547240956, 1.297916304, 1.111509602, 0.951874626, 
    0.806784197, 0.656114181, 0.55038699, 0.46649371, 0.420684793, 0.383315809, 0.367791036, 
    0.364009101, 0.367791036, 0.379374236, 0.391322237, 0.395387946, 0.399495897, 0.383315809, 
    0.342120308, 0.30221226, 0.256147259, 0.21710376, 0.187855007, 0.171168046, 0.162546717, 
    0.162546717, 0.169407953, 0.180246642, 0.191778788, 0.204048759, 0.214871318, 0.226267895, 
    0.235818856, 0.245772972, 0.256147259, 0.269733074, 0.284039469, 0.296029017, 0.305352146,
    0.311730124, 0.318241321, 0.318241321, 0.318241321};
    
    double m_CD_sphere = 0.41;

    // declare variables for arrow_rk4
    double* m_k1 = new double[12];
    double* m_k2 = new double[12];
    double* m_k3 = new double[12];
    double* m_k4 = new double[12];
    double* m_y_temp = new double[12];


    // functions;

    void aerodynamics_cannonball(double* y, double* ans);
    double get_sphere_CD(double reynolds);

    void cannonball_rk4_func(double t, double* y, double* ans);
    void cannonball_rk4(double t0, double* y0, double dt, int size, double* ans);

    

};

#endif