//#include "hoch_helper.h"
#include "aircraft.h"

int main (int argc, char * const argv[])
{
    cout<<"Program started: \n"<<endl;

    
    double  euler[3];
    double  quat[4];
    euler[0] = 0.0;
    euler[1] = pi/2.0;
    euler[2] = 0.0;
    euler_to_quat(&euler[0], &quat[0]);
    array_print(quat, 4);
    
    aircraft* my_aircraft = new aircraft("vertical_trim2.json");
    my_aircraft->init_sim();
    printf("Running sim\n\n");
    my_aircraft->run_sim();

    
    cout<<"Program completed: "<<endl;

    return 0;
}