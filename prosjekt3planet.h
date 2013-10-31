#ifndef PLANET_H
#define PLANET_H
#include <armadillo>


using namespace std;
using namespace arma;

class Planet
{

public:
    Planet(string inname, double inmass);

    Planet();

    string name;
    double mass;

    double xpos;
    double ypos;
    double xvel;
    double yvel;

    void setinitialconditions(double xpos, double ypos, double xvel, double yvel);
};

#endif // PLANET_H
