
#include "planet.h"
#include "solsystem.h"
#include <math.h>
#include <armadillo>
#include <iostream>

using namespace arma;

//only used to find initial velocity of sun
double momentum(Planet* inplanets,int n){
    double momentum = 0;
    for(int i = 0; i<=n-1; i++){
        momentum += inplanets[i].mass*inplanets[i].xvel;
    }
    return momentum;
}

int main()
{
    double refmass = 2*pow(10.0,30.0);
    double refpos = 149.6; //10^6 KM | all planets on a line at start.
    double refvel = 29.8/(2*acos(-1)); // earthvelocity/2pi should give one

    //mercury
    double mercurymass = 0.33*pow(10.0,24.0)/refmass;
    Planet* mercury = new Planet("mercury", mercurymass);
    double mercurypos = 57.9/refpos;
    double mercuryvel = 47.9/refvel;
    mercury->setinitialconditions(0,mercurypos,mercuryvel,0);

    //venus
    double venusmass = 4.87*pow(10.0,24.0)/refmass;
    Planet* venus = new Planet("venus", venusmass);
    double venuspos = 108.2/refpos;
    double venusvel = 35.0/refvel;
    venus->setinitialconditions(0,venuspos,venusvel,0);

    //earth
    double earthmass = 6*pow(10.0,24.0)/refmass;
    Planet* earth = new Planet("earth", earthmass);
    double earthvel = 2*acos(-1);
    earth->setinitialconditions(0,1,earthvel,0);

    //mars
    double marsmass = 0.642*pow(10.0,24.0)/refmass;
    Planet* mars = new Planet("mars", marsmass);
    double marspos = 227.9/refpos;
    double marsvel = 24.1/refvel;
    mars->setinitialconditions(0,marspos,marsvel,0);

    //jupiter
    double jupitermass = 1.9*pow(10.0,24.0)/refmass;
    Planet* jupiter = new Planet("jupiter", jupitermass);
    double jupiterpos = 778.6/refpos;
    double jupitervel = 13.1/refvel;
    jupiter->setinitialconditions(0,jupiterpos,jupitervel,0);

    //saturn
    double saturnmass = 568*pow(10.0,24.0)/refmass;
    Planet* saturn = new Planet("saturn", saturnmass);
    double saturnpos = 1433.5/refpos;
    double saturnvel = 9.7/refvel;
    saturn->setinitialconditions(0,saturnpos,saturnvel,0);

    //uranus
    double uranusmass = 86.8*pow(10.0,24.0)/refmass;
    Planet* uranus = new Planet("uranus", uranusmass);
    double uranuspos = 2872.5/refpos;
    double uranusvel = 6.8/refvel;
    uranus->setinitialconditions(0,uranuspos,uranusvel,0);

    //neptune
    double neptunemass = 102*pow(10.0,24.0)/refmass;
    Planet* neptune = new Planet("neptune", neptunemass);
    double neptunepos = 4495.1/refpos;
    double neptunevel = 5.4/refvel;
    neptune->setinitialconditions(0,neptunepos,neptunevel,0);

    //pluto
    double plutomass = 0.0131*pow(10.0,24.0)/refmass;
    Planet* pluto = new Planet("pluto", plutomass);
    double plutopos = 5870.0/refpos;
    double plutovel = 4.7/refvel;
    pluto->setinitialconditions(0,plutopos,plutovel,0);

    //sun
    double sunmass = 1;
    Planet* sun = new Planet("sun", sunmass);
    sun->setinitialconditions(0,0,0,0);//set to zero to find momentum with momentummethod


    int input = 0;

    cout << "1) earth and sun" << endl;
    cout << "2) earth, sun and jupiter" << endl;
    cout << "3) from sun to mars" << endl;
    cout << "4) all planets" << endl;

    cin >> input;
    Planet* inplanets;
    double sunvel; //velocity of sun to get zero momentum
    int n; //number of planets
    switch(input){
        case 1:
            n = 2;
            inplanets = new Planet[n];
            inplanets[0] = *sun;
            inplanets[1] = *earth;
            sunvel = -momentum(inplanets,n)/sunmass;
            sun->setinitialconditions(0,0,sunvel,0);
            inplanets[0] = *sun;

            break;
        case 2:
            n = 3;
            inplanets = new Planet[n];
            inplanets[0] = *sun;
            inplanets[1] = *earth;
            inplanets[2] = *jupiter;
            sunvel = -momentum(inplanets,n)/sunmass;
            sun->setinitialconditions(0,0,sunvel,0);
            inplanets[0] = *sun;
            break;
        case 3:
            n = 5;
            inplanets = new Planet[n];
            inplanets[0] = *sun;
            inplanets[1] = *mercury;
            inplanets[2] = *venus;
            inplanets[3] = *earth;
            inplanets[4] = *mars;
            sunvel = -momentum(inplanets,n)/sunmass;
            sun->setinitialconditions(0,0,sunvel,0);
            inplanets[0] = *sun;
            break;
        case 4:
            n = 10;
            inplanets = new Planet[n];
            inplanets[0] = *sun;
            inplanets[1] = *mercury;
            inplanets[2] = *venus;
            inplanets[3] = *earth;
            inplanets[4] = *mars;
            inplanets[5] = *jupiter;
            inplanets[6] = *saturn;
            inplanets[7] = *neptune;
            inplanets[8] = *uranus;
            inplanets[9] = *pluto;
            sunvel = -momentum(inplanets,n)/sunmass;
            sun->setinitialconditions(0,0,sunvel,0);
            inplanets[0] = *sun;
            break;
    }
    Solsystem melkeveien(inplanets,n);

    double time = 0;
    cout << "time in years" << endl;
    cin >> time;

    double h = 1;
    cout << "timestep" << endl;
    cin >> h;

    delete[] inplanets;

    cout << "momentum:" << momentum(inplanets,n) << endl;

    cout << "before:" << melkeveien.potensialenergy() + melkeveien.kineticenergy() << endl;

    melkeveien.movement(time,h);
    cout << "after:" << melkeveien.potensialenergy() + melkeveien.kineticenergy() << endl;

    return 0;
}
