//#include "hoch_helper.h"
#include "aircraft.h"

int main (int argc, char * const argv[])
{
    cout<<"Program started: \n"<<endl;

    aircraft* my_aircraft = new aircraft("aircraft_7_7_4.json");
    
    my_aircraft->init_sim();
    printf("Running sim\n\n");
    my_aircraft->run_sim();

    
    cout<<"Program completed: "<<endl;

    return 0;
}