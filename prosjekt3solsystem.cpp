#include "solsystem.h"
#include "planet.h"
#include <math.h>
#include <stdio.h>


Solsystem::Solsystem(Planet* inplanets, int innumberofplanets){
    fourpi2 = 4*acos(-1)*acos(-1);
    numberofplanets = innumberofplanets;
    planets = new Planet[numberofplanets];
    for(int i = 0; i<= numberofplanets-1; i++){
        planets[i] = inplanets[i];
    }
    initialenergy = kineticenergy() + potensialenergy();
}

void Solsystem::printoutplanets(void){
    for (int i = 0; i<= numberofplanets - 1; i++){
        cout << planets[i].name << endl;
    }

}

void Solsystem::printoutpositions(void){
    for(int i = 0; i<= numberofplanets - 1; i++){
        cout << planets[i].name << endl;
        cout << planets[i].xpos << "   " << planets[i].ypos << endl << endl;
        cout << sqrt(pow(planets[i].xpos,2.0) + pow(planets[i].ypos, 2.0)) << endl << endl;
    }
}

double Solsystem::potensialenergy(void){
    double energy = 0;
    double G = 4*acos(-1)*acos(-1);
    double r;
    for (int i = 0; i<=numberofplanets-1; i++){
        for (int j = i+1; j<=numberofplanets-1; j++){
            r = pow(planets[i].xpos-planets[j].xpos,2.0) + pow(planets[i].ypos-planets[j].ypos,2.0);
            energy += G*planets[j].mass*planets[i].mass/r;
        }
    }
    return energy;
}

double Solsystem::kineticenergy(void){
    double energy = 0;
    double v;
    for(int i = 0; i<=numberofplanets-1; i++){
        v = pow(planets[i].xvel,2) + pow(planets[i].yvel,2);
        energy += 0.5*planets[i].mass*v*v;
    }
    return energy;
}

void Solsystem::printtofile(void){
    ofstream fil;
    fil.open("positions.txt", ofstream::app);
    fil << fixed;
    fil.precision(5);
    for (int i = 0; i<=numberofplanets-1; i++){
        fil << planets[i].xpos << " " << planets[i].ypos << " ";
    }
    fil << endl;
    fil.close();

    fil.open("energy.txt", ofstream::app);
    fil << fixed;
    fil.precision(5);
    double energydiff;
    for (int i = 0; i<=numberofplanets-1; i++){
        energydiff = initialenergy - (kineticenergy() + potensialenergy());
        fil << energydiff << endl;
    }
    fil.close();
}

// force for one planet in booth coordinates
void Solsystem::force(int planet, double* xforce, double* yforce){
    *xforce = 0;
    *yforce = 0;
    double r;
    double planetxpos = planets[planet].xpos;
    double planetypos = planets[planet].ypos;
    for (int i = 0; i<=numberofplanets-1; i++){
        if(i!=planet){
            r = sqrt(pow(planets[i].xpos-planetxpos,2.0) + pow(planets[i].ypos-planetypos,2.0));
            *xforce += fourpi2*planets[i].mass*((planets[i].xpos-planetxpos)/(r*r*r));
            *yforce += fourpi2*planets[i].mass*((planets[i].ypos-planetypos)/(r*r*r));
        }
    }
}

void Solsystem::movement(double intime, double inh){
    double time = intime;
    double h = inh;
    double step = time/h;
    remove("positions.txt");
    remove("energy.txt");
    for(int i = 0; i<= step; i++){
        RK4(h);
        if(i%10==0){
            printtofile();
        }
        /*if(i%1000000==0){
            cout << i << endl;
        }*/
    }
}

void Solsystem::RK4(double h){
    //array for storing inintial velocity
    double *initxvel = new double[numberofplanets];
    double *inityvel = new double[numberofplanets];

    //array for storing initial posistion
    double *initxpos = new double[numberofplanets];
    double *initypos = new double[numberofplanets];

    for (int i = 0; i<= numberofplanets-1; i++){
        initxvel[i] = planets[i].xvel;
        inityvel[i] = planets[i].yvel;
        initxpos[i] = planets[i].xpos;
        initypos[i] = planets[i].ypos;
    }


    //array for stroing of k
    double *xvelk = new double[4*numberofplanets];
    double *yvelk = new double[4*numberofplanets];
    double *xposk = new double[4*numberofplanets];
    double *yposk = new double[4*numberofplanets];

    for (int i = 0; i<=numberofplanets-1; i++){
        double forcex = 0;
        double forcey = 0;

        //0 (k1s)
        force(i, &forcex, &forcey);

        xvelk[4*i] = h*forcex;
        yvelk[4*i] = h*forcey;

        xposk[4*i] = h*initxvel[i];
        yposk[4*i] = h*inityvel[i];
        //end 0 (k1s)

        //1 (k2s)
        //eulers formula to find halfsteps
        planets[i].xvel = initxvel[i] + xvelk[4*i]/2.0;
        planets[i].yvel = inityvel[i] + yvelk[4*i]/2.0;

        planets[i].xpos = initxpos[i] + xposk[4*i]/2.0;
        planets[i].ypos = initypos[i] + yposk[4*i]/2.0;

        //new forces for halfstep
        force(i, &forcex, &forcey);

        xvelk[4*i+1] = h*forcex;
        yvelk[4*i+1] = h*forcey;

        xposk[4*i+1] = h*(initxvel[i]+xvelk[4*i]/2.0);
        yposk[4*i+1] = h*(inityvel[i]+yvelk[4*i]/2.0);
        //end 1 (k2s)

        //2 (k3s)
        //eulers formula to find halfsteps
        planets[i].xvel = initxvel[i] + xvelk[4*i+1]/2.0;
        planets[i].yvel = inityvel[i] + yvelk[4*i+1]/2.0;

        planets[i].xpos = initxpos[i] + xposk[4*i+1]/2.0;
        planets[i].ypos = initypos[i] + yposk[4*i+1]/2.0;

        //new forces for halfstep
        force(i, &forcex, &forcey);

        xvelk[4*i+2] = h*forcex;
        yvelk[4*i+2] = h*forcey;

        xposk[4*i+2] = h*(initxvel[i]+xvelk[4*i+1]/2.0);
        yposk[4*i+2] = h*(inityvel[i]+yvelk[4*i+1]/2.0);

        //end 2 (k3s)

        //3 (k4s)
        //eulers formula to find halfsteps
        planets[i].xvel = initxvel[i] + xvelk[4*i+2];
        planets[i].yvel = inityvel[i] + yvelk[4*i+2];

        planets[i].xpos = initxpos[i] + xposk[4*i+2];
        planets[i].ypos = initypos[i] + yposk[4*i+2];

        //new forces for halfstep
        force(i, &forcex, &forcey);

        xvelk[4*i+3] = h*forcex;
        yvelk[4*i+3] = h*forcey;

        xposk[4*i+3] = h*(initxvel[i] + xvelk[4*i+2]);
        yposk[4*i+3] = h*(inityvel[i] + yvelk[4*i+2]);

        //end 3 (k4s)
    }

    for(int i = 0; i<=numberofplanets-1; i++){
        planets[i].xvel = initxvel[i] + (xvelk[4*i+0] + 2*xvelk[4*i+1] + 2*xvelk[4*i+2] + xvelk[4*i+3])/6.0;
        planets[i].yvel = inityvel[i] + (yvelk[4*i+0] + 2*yvelk[4*i+1] + 2*yvelk[4*i+2] + yvelk[4*i+3])/6.0;

        planets[i].xpos = initxpos[i] + (xposk[4*i+0] + 2*xposk[4*i+1] + 2*xposk[4*i+2] + xposk[4*i+3])/6.0;
        planets[i].ypos = initypos[i] + (yposk[4*i+0] + 2*yposk[4*i+1] + 2*yposk[4*i+2] + yposk[4*i+3])/6.0;
    }
}
