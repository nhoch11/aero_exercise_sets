//#include "hoch_helper.h"
#include "aircraft.h"

int main (int argc, char * const argv[])
{
    cout<<"Program started: "<<endl;

    aircraft* my_aircraft = new aircraft("aircraft_6_4.json");

    my_aircraft->init_sim();
    my_aircraft->exercise_6_4();

    
    cout<<"Program completed: "<<endl;

    return 0;
}