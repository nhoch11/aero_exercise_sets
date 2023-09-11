#include "hoch_helper.h"
#include "arrow.h"

int main (int argc, char * const argv[])
{
    arrow* my_arrow = new arrow("arrow.json");

    my_arrow->shoot_arrow();

    return 0;
}