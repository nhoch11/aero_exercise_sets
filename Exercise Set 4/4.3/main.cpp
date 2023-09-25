//#include "hoch_helper.h"
#include "arrow.h"

int main (int argc, char * const argv[])
{
    cout<<"Program started: "<<endl;

    arrow* my_arrow = new arrow("arrow.json");

    my_arrow->exercise_3_9();


    cout<<"Program completed: "<<endl;

    return 0;
}