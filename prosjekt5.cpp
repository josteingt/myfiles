#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

void eksplisitt(int ntime, int npos, double alpha, double* unew, double* uold){

    for(int i = 1; i <= ntime-2; i++){
        for( int j = 1; j<= npos-2; j++){
            unew[j] = uold[j] + alpha*(uold[j+1] + uold[j-1] - 2*uold[j]);
        }
        for(int k = 0; k<=npos-1; k++){
            uold[k] = unew[k];
        }
    }

}

void implisitt(int ntime, int npos, double alpha, double* unew, double* uold){
    double* diagonal = new double[npos];
    double* nondiagonalupper = new double[npos-1];
    double* nondiagonallower = new double[npos-1];
    double diagonalterm = 1 + 2*alpha;

    for(int j = 1; j <= ntime - 1; j++){

        //setting up matrix
        for(int i = 0; i<= npos-2; i++){
            diagonal[i] = diagonalterm;
            nondiagonallower[i] = -alpha;
            nondiagonalupper[i] = -alpha;
        }

        diagonal[npos-1] = diagonalterm;
        uold[1] += alpha;
        //forward substitution
        double temp = 0;
        for(int i = 2; i<=npos-2;i++){
            temp = nondiagonallower[i]/diagonal[i-1];
            diagonal[i] -= temp*nondiagonalupper[i];
            uold[i] -= temp*uold[i-1];
            nondiagonallower[i] = 0;
        }

//        cout << endl;
//        for (int i = 0; i<= npos -1; i++){
//            cout << "i=" << i << "   ";
//            cout << nondiagonallower[i] << "   ";
//            cout << diagonal[i] << "   ";
//            cout << nondiagonalupper[i] << "   ";
//            cout <<  uold[i] << endl;
//        }
//        cout << endl;

        //backward substitusjon
        for (int i = npos-2; i>=1; i--){
            uold[i] = (uold[i] - nondiagonalupper[i+1]*uold[i+1])/diagonal[i];
            nondiagonalupper[i+1] = 0;
            diagonal[i] = 1;
        }
    }

}

void cranknicolson(int ntime, int npos, double alpha, double* unew, double* uold){

    double* diagonal = new double[npos];
    double* nondiagonalupper = new double[npos];
    double* nondiagonallower = new double[npos];
    double diagonalterm = 2 + 2*alpha;

    for(int j = 1; j <= ntime - 1; j++){

        //setting up matrix
        for(int i = 0; i<= npos-1; i++){
            diagonal[i] = diagonalterm;
            nondiagonallower[i] = -alpha;
            nondiagonalupper[i] = -alpha;
        }
        diagonal[npos-1] = diagonalterm;

        //multiplication with 2I-alphaB
        for(int i = 1; i<= npos-2; i++){
            unew[i] = alpha*uold[i-1] + (2-2*alpha)*uold[i] + alpha*uold[i+1];
        }
        for(int i = 0; i<= npos-1; i++){
            uold[i] = unew[i];
        }

        uold[1] += alpha; //because of initial conditions

        //forward substitution
        double temp = 0;
        for(int i = 2; i<=npos-2;i++){
            temp = nondiagonallower[i]/diagonal[i-1];
            diagonal[i] -= temp*nondiagonalupper[i];
            uold[i] -= temp*uold[i-1];
            nondiagonallower[i] = 0;
        }

//        cout << endl;
//        for (int i = 0; i<= npos -1; i++){
//            cout << "i=" << i << "   ";
//            cout << nondiagonallower[i] << "   ";
//            cout << diagonal[i] << "   ";
//            cout << nondiagonalupper[i] << "   ";
//            cout <<  uold[i] << endl;
//        }
//        cout << endl;

        //backward substitusjon
        for (int i = npos-2; i>=1; i--){
            uold[i] = (uold[i] - nondiagonalupper[i+1]*uold[i+1])/diagonal[i];
            nondiagonalupper[i+1] = 0;
            diagonal[i] = 1;
        }
    }
}

int main()
{
    //setting position
    double a = 0;
    double b = 1;

    int npos = 10;
    double deltax = (b-a)/(npos-1);

    double alpha = 0.5;

    //calculating time
    double time = 10;
    double deltat = alpha*deltax*deltax;
    int ntime = time/deltat + 2;// +1 for correct amount +1 to be safe with truncation by int.

    cout << "ntime:  " << ntime << endl;
    cout << "npos:   " << npos << endl;

    //explicit scheme
    //setting up inintial condintions
    double* unew = new double[npos];
    double* uold = new double[npos];

    for (int i = 0; i<= npos-1; i++){
        unew[i] = 0;
        uold[i] = 0;
    }

    unew[0] = 1;
    uold[0] = 1;

    cout << "initial conditions:" << endl;
    for (int i = 0; i<= npos -1; i++){
        cout << "i=" << i << "   " << uold[i] << endl;
    }
    cout << endl;
    cout << endl << "Explicit:" << endl;

    //solving
    eksplisitt(ntime, npos, alpha, unew, uold);

    cout << "solution:" << endl;
    for (int i = 0; i<= npos -1; i++){
        cout << "i=" << i << "   " << uold[i] << endl;
    }
    cout << endl;
    //end explicit scheme

    ofstream fil;
    fil.open("eksplisitt.txt");
    fil << fixed;
    fil.precision(5);
    for(int i = 0; i<=npos-1; i++){
        fil << i*deltax << " " << uold[i] << endl;
    }
    fil.close();

    //implicit scheme
    //setting up inintial condintions
    for (int i = 0; i<= npos-1; i++){
        unew[i] = 0;
        uold[i] = 0;
    }

    unew[0] = 1;
    uold[0] = 1;

    cout << endl << "Implicit:" << endl;
    //solving
    implisitt(ntime,npos,alpha,unew,uold);

    cout << "solution:" << endl;
    for (int i = 0; i<= npos -1; i++){
        cout << "i=" << i << "   " << uold[i] << endl;
    }
    cout << endl;
    //end implicit scheme

    fil.open("implisitt.txt");
    fil << fixed;
    fil.precision(5);
    for(int i = 0; i<=npos-1; i++){
        fil << i*deltax << " " << uold[i] << endl;
    }
    fil.close();

    //crank-nicolson scheme
    //setting up inintial condintions
    for (int i = 0; i<= npos-1; i++){
        unew[i] = 0;
        uold[i] = 0;
    }

    unew[0] = 1;
    uold[0] = 1;

    cout << endl << "Crank nicolson:" << endl;

    //solving
    cranknicolson(ntime,npos,alpha,unew,uold);

    cout << "solution:" << endl;
    for (int i = 0; i<= npos -1; i++){
        cout << "i=" << i << "   " << uold[i] << endl;
    }
    cout << endl;
    //end crank-nicolson

    fil.open("cn.txt");
    fil << fixed;
    fil.precision(5);
    for(int i = 0; i<=npos-1; i++){
        fil << i*deltax << " " << uold[i] << endl;
    }
    fil.close();
    //Analytical

    double* analytical = new double[npos];

    for(int i = 0; i<= npos-1; i++){
        analytical[i] = 0;
    }

    int numberofterms = 10;
    double pi = acos(-1);
    for(int i = 0; i <= npos-1; i++){// position
        for(int j = 1; j <= numberofterms; j++){// sum of sin()
            analytical[i] -= 2/(j*pi)*sin(j*pi*i*deltax)*exp(-j*j*pi*pi*time);
        }
        analytical[i] += 1-i*deltax;
    }
    cout << "solution:" << endl;
    for (int i = 0; i<= npos -1; i++){
        cout << "i=" << i << "   " << analytical[i] << endl;
    }

    fil.open("analytical.txt");
    fil << fixed;
    fil.precision(5);
    for(int i = 0; i<=npos-1; i++){
        fil << i*deltax << " " << analytical[i] << endl;
    }
    fil.close();
}

