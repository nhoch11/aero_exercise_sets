//#include "hoch_helper.h"
#include "cannonball.h"

int main (int argc, char * const argv[])
{
    cout<<"Program started: "<<endl;

    cannonball* my_cannonball = new cannonball("cannonball.json");

    my_cannonball->exercise_4_3();

    
    cout<<"Program completed: "<<endl;

    return 0;
}