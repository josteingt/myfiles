#ifndef SOLSYSTEM_H
#define SOLSYSTEM_H



#include "planet.h"

class Solsystem{
private:
    double initialenergy;

    int numberofplanets;

    double fourpi2;

    Planet *planets;

    void RK4(double h);

    void force(int, double*, double*);

    void printtofile(void);

public:

    Solsystem(Planet* inplanets, int numberofplanets);

    void printoutplanets(void);

    void printoutpositions(void);

    void movement(double time, double h);

    double potensialenergy(void);

    double kineticenergy(void);

};

#endif // SOLSYSTEM_H
