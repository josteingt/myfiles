#include "planet.h"
#include <math.h>
#include <armadillo>

using namespace arma;
using namespace std;



//constructor name and mass
Planet::Planet(string inname, double inmass){
    name = inname;
    mass = inmass;
}

Planet::Planet(){

}


//setting initialconditions for position and velocity
void Planet::setinitialconditions(double inxpos, double inypos, double inxvel, double inyvel){
    xpos = inxpos;
    ypos = inypos;

    xvel = inxvel;
    yvel = inyvel;
}
