//#include "hoch_helper.h"
#include "aircraft.h"

int main (int argc, char * const argv[])
{
    cout<<"Program started: \n"<<endl;

    //aircraft* my_aircraft = new aircraft("aircraft_6_6_test.json");
    
    //my_aircraft->init_sim();
    //printf("Running sim\n\n");
    //my_aircraft->run_sim();
    double A[4], B[4], ans[4];
    A[0] = 1.0;
    A[1] = 2.0;
    A[2] = 3.0;
    A[3] = 4.0;
    B[0] = 1.0;
    B[1] = -2.0;
    B[2] = -3.0;
    B[3] = -4.0;
    quat_mult(A, B, ans);
    array_print(ans, 4);

    cout<<"Program completed: "<<endl;

    return 0;
}