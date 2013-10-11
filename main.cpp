#include <iostream>
#include <armadillo>
#include <math.h>
#include <time.h>

using namespace std;
using namespace arma;

double absverdi(double test){
    double retur;
    if(test>=0){
        retur = test;
    }
    else{
        retur = -test;
    }
    return retur;
}

double offdiag(mat& matrise, int* maxrow, int* maxcol, int size){
    double max = 0;
    double aij;

    for (int i = 0; i<=size-1; i++){
        for(int j = i + 1; j<=size-1; j++){
            aij = absverdi(matrise(i,j));
            if(aij > max){
                max = aij;
                *maxrow = i;
                *maxcol = j;
            }
        }
    }
    return max;
}

void jacobirotate(mat& schrodinger, mat& egenvektor, int maxrow, int maxcol, int size){
    double sin, cos;
    if(schrodinger(maxcol,maxrow) !=0){
        double t, tau;
        tau = (schrodinger(maxrow,maxrow) - schrodinger(maxcol,maxcol))/(2*schrodinger(maxcol,maxrow));
        if(tau>0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else{
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }
        cos = 1.0/sqrt(1.0+t*t);
        sin = cos*t;
    }
    else{
        cos = 1;
        sin = 0;
    }

    double a_colcol, a_rowrow, a_icol, a_irow, r_icol, r_irow;
    a_colcol = schrodinger(maxcol,maxcol);
    a_rowrow = schrodinger(maxrow,maxrow);

    //changing elements with maxcol and maxrow
    schrodinger(maxcol,maxcol) = cos*cos*a_colcol - 2.0*cos*sin*schrodinger(maxcol,maxrow) + sin*sin*a_rowrow;
    schrodinger(maxrow,maxrow) = sin*sin*a_colcol + 2.0*cos*sin*schrodinger(maxcol,maxrow) + cos*cos*a_rowrow;
    schrodinger(maxcol,maxrow) = 0;
    schrodinger(maxrow,maxcol) = 0;

    for(int i = 0; i<=size-1; i++){
        if(i != maxcol && i != maxrow){
            a_icol = schrodinger(i,maxcol);
            a_irow = schrodinger(i,maxrow);
            schrodinger(i,maxcol) = cos*a_icol - sin*a_irow;
            schrodinger(maxcol,i) = schrodinger(i,maxcol);
            schrodinger(i,maxrow) = cos*a_irow + sin*a_icol;
            schrodinger(maxrow,i) = schrodinger(i,maxrow);
        }
        //changing eigenvector
        r_icol = egenvektor(i,maxcol);
        r_irow = egenvektor(i,maxrow);
        egenvektor(i,maxcol) = cos*r_icol - sin*r_irow;
        egenvektor(i,maxrow) = cos*r_irow + sin*r_icol;
    }
}

void jacobimetode(mat& schrodinger, mat& egenvektor, int size){
    int maxrow, maxcol;
    double precision = pow(10.0,-8.0);
    int maxiter = size*size*size;
    int iterations = 0;

    double maxoffdiag = offdiag(schrodinger,&maxrow,&maxcol,size);

    while(maxoffdiag > precision && iterations < maxiter){
        maxoffdiag = offdiag(schrodinger, &maxrow, &maxcol, size);
        jacobirotate(schrodinger, egenvektor, maxrow, maxcol, size);
        iterations++;
        //cout << schrodinger << endl;
    }
    cout << "iterations: " << endl << iterations << endl;

}

int main(){

    double pmax = 5.0;
    double pmin = 0.0;
    int steg = 202;
    double omega = 5;
    int matstor = steg - 2;
    double h = (pmax-pmin)/(steg-1);
    double hh = h*h;

    vec p = vec(matstor);
    for (int i = 0; i<=matstor-1; i++){
        p(i) = pmin + (i+1)*h;
    }

    //inintialize schrodinger
    mat schrodinger = mat(matstor,matstor);
    for(int i = 0; i<=matstor-1; i++){
        for(int j = 0; j<= matstor-1; j++){
            if(i == j){
                schrodinger(i,j) = 2.0/hh + omega*omega*p(j)*p(j);
            }
            else if(i == j-1 || i-1 == j){
                schrodinger(i,j) = -1.0/hh;
            }
            else{
                schrodinger(i,j) = 0.0;
            }
        }
    }
    mat copyschrodinger = schrodinger;
    //cout << schrodinger << endl;

    // initialize eigenvector
    mat egenvektor = eye<mat>(matstor,matstor);


    // running jacobimethod
    clock_t t = clock();
    jacobimetode(schrodinger, egenvektor, matstor);
    cout << "jacobi solver: " << ((double)(clock() - t))/CLOCKS_PER_SEC << endl;

    // test if vectors are normalized
    /*double normfaktor;
    for (int i = 0; i <=matstor-1; i++){
        normfaktor = 0;
        for(int j = 0; j<=matstor-1; j++){
            normfaktor += egenvektor(i,j)*egenvektor(i,j);
        }
        cout << normfaktor << endl;
    }*/
    // vectors are normalized out of the algorithm

    //finding index of lowest eigenvalue
    double min = 10000;
    int minindex;
    for(int i = 0; i<= matstor-1; i++){
        if(schrodinger(i,i) < min){
            min = schrodinger(i,i);
            minindex = i;
        }
    }

    //normalizing for the integral
    //creating psi^2
    vec psikvadrat = vec(matstor);
    for(int i = 0; i<=matstor-1; i++){
        psikvadrat(i) = egenvektor(i,minindex)*egenvektor(i,minindex);
    }

    //normalizing
    for (int i = 0; i<= matstor-1; i++){
        psikvadrat(i) /= h;
    }

    //cout << egenvektor << endl;
    //extracting eigenvalues
    vec eigschrodinger = sort(schrodinger.diag());

    cout << "mine:" << endl << eigschrodinger(0) << endl << eigschrodinger(1) << endl << eigschrodinger(2) << endl;
    //cout << egenvektor << endl;

    cx_vec eigval;
    cx_mat eigvec;
    t = clock();
    eig_gen(eigval,eigvec,copyschrodinger);
    cout << "armadillo: " << ((double)(clock() - t))/CLOCKS_PER_SEC << endl;
    eigval = sort(eigval);
    cout << "armadillo sine: " << endl << eigval(0) << endl << eigval(1) << endl << eigval(2) << endl;
    //cout << eigvec << endl;

    //write psikvadrat to file

    ofstream fil;
    fil.open("psikvadrat5without.txt");
    fil << fixed;
    fil.precision(5);
    fil << 0 << " " << 0 << endl;
    for(int i = 0; i<matstor-1; i++){
        fil << p(i) << " " << psikvadrat(i) << endl;
    }
    fil.close();

}
