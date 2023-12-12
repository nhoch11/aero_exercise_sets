//#include "hoch_helper.h"
#include "aircraft.h"

int main (int argc, char * const argv[])
{
    cout<<"Program started: \n"<<endl;

<<<<<<< HEAD
    aircraft* my_aircraft = new aircraft("aircraft_7_7_3.json");
=======
    aircraft* my_aircraft = new aircraft("aircraft_6_6_test.json");
>>>>>>> c375a94b9fab333801684dcbabc26d7154b2fd94
    
    my_aircraft->init_sim();
    printf("Running sim\n\n");
    my_aircraft->run_sim();

    
    cout<<"Program completed: "<<endl;

    return 0;
}