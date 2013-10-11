#include <iostream>
#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

double absverdi(double test){
    double retur;
    if(test >= 0){
        retur = test;
    }
    else{
        retur = -test;
    }
    return retur;
}

double offdiag(mat* matrise, int* maxrow, int* maxcol, int size){
    double max = 0;
    double aij;
    for(int i = 0; i<=size-1; i++){
        for(int j = i + 1; j<=size-1; j++){
            if(i != j){
                aij = absverdi((*matrise)(i,j));
                if(aij>max){
                    max = aij;
                    *maxrow = i;
                    *maxcol = j;
                }
            }
        }
    }
    cout << *maxrow << endl << *maxcol << endl;
    return max;
}

void jacobi(mat* egenverdi, mat* egenvektor,int maxrow, int maxcol, int size){
    double sin, cos;
    if((*egenverdi)(maxrow,maxcol) != 0){
        double tau = ((*egenverdi)(maxcol,maxcol) - (*egenverdi)(maxrow,maxrow))/2*(*egenverdi )(maxrow,maxcol);

        double t;
        if (tau>0){
            t = 1.0/(tau + sqrt(1+tau*tau));
        }
        else{
            t = -1.0/(-tau + sqrt(1+tau*tau));
        }
        cos = 1/sqrt(1+t*t);
        sin = t*cos;
    }
    else{
        cos = 1.0;
        sin = 0.0;
    }
    double a_rowrow, a_colcol, a_irow, a_icol, r_irow, r_icol;
    a_rowrow= (*egenverdi)(maxrow,maxrow);
    a_colcol = (*egenverdi)(maxcol,maxcol);

    // changigng elements with indices maxrow and maxcol
    (*egenverdi)(maxrow,maxrow) = cos*cos*a_rowrow - 2.0*cos*sin*(*egenverdi)(maxrow,maxcol) + sin*sin*a_colcol;
    (*egenverdi)(maxcol,maxcol) = sin*sin*a_rowrow + 2.0*cos*sin*(*egenverdi)(maxrow,maxcol) + cos*cos*a_colcol;
    (*egenverdi)(maxrow,maxcol) = 0;
    (*egenverdi)(maxcol,maxrow) = 0;

    // canging the rest
    for(int i = 0; i<=size-1; i++){
        if(i != maxrow && i != maxcol){
            a_irow = (*egenverdi)(i,maxrow);
            a_icol = (*egenverdi)(i,maxcol);
            (*egenverdi)(i,maxrow) = cos*a_irow - sin*a_icol;
            (*egenverdi)(maxrow,i) = (*egenverdi)(i,maxrow);
            (*egenverdi)(i,maxcol) = cos*a_icol + sin*a_irow;
            (*egenverdi)(maxcol,i) = (*egenverdi)(i,maxcol);
        }
        //changing the eigenvectors
        r_irow = (*egenvektor)(i,maxrow);
        r_icol = (*egenvektor)(i,maxcol);
        (*egenvektor)(i,maxrow) = cos*r_irow - sin*r_icol;
        (*egenvektor)(i,maxcol) = sin*r_irow + cos*r_icol;
    }
    return;
}

int main()
{
    double pmax = 10.0;
    double pmin = 0.0;
    int steg = 20;
    double h = (pmax - pmin)/(steg-1);
    double hh = h*h;

    //Posistion vector
    mat p = mat(1,steg);
    for(int i = 0; i<=steg-1; i++){
        p(i) = pmin + i*h;
    }

    // initialize matrix
    mat schrodinger = mat(steg,steg);

    for(int i = 0; i<=steg-1; i++){
        for(int j =  0; j<= steg-1; j++){
            if ( i == j){
                schrodinger(i,j) = 2/(hh) + p(i)*p(i);
            }
            else if(i == j-1 || i-1 == j){
                schrodinger(i,j) = -1/hh;
            }
            else{
                schrodinger(i,j) = 0;
            }
        }
    }
    //inintialize eigentvectormatrix
    mat egenvektor = eye<mat>(steg,steg);


    mat copyschrodinger = schrodinger;
    double maxnondiag = 10;
    int maxiter = 5;
    double tolerance = pow(10.0,-10.0);
    int iterations = 0;
    while (maxnondiag > tolerance && iterations < maxiter){
        int maxrow, maxcol;
        maxnondiag = offdiag(&schrodinger, &maxrow, &maxcol, steg);
        jacobi(&schrodinger,&egenvektor,maxrow,maxcol,steg);
        //cout << schrodinger << endl;
        iterations++;
    }

    //egenverdier i vektor
    vec eigschrodinger = vec(steg);
    for (int i = 0; i<=steg-1; i++){
        eigschrodinger(i) = schrodinger(i,i);
    }
    cout << " mine:" << endl << eigschrodinger << endl;
    //cout << iterations << endl;

    cx_vec eigval;
    cx_mat eigvec;
    eig_gen(eigval,eigvec,copyschrodinger);
    cout << "armadillo sine: " << endl << eigval << endl;
    return 0;
}

