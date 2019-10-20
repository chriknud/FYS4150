#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include "armadillo"
#include "jacobi.h"
# include <assert.h>
using namespace std;
using namespace arma;

//  we have defined a matrix A and a matrix R for the eigenvector, both of dim n x n
//  The final matrix R has the eigenvectors in its row elements, it is set to one
//  for the diagonal elements in the beginning, zero else.

int main()
{
    jacobi jacobi_method;
    int n = 200;
    mat A(n,n,fill::zeros);
    mat R(n,n,fill::eye);
    double tolerance = 1.0E-15;
    int iterations = 0;
    int maxiter = 1000000;
    double maxnondiag = tolerance*2;
    double rho0 = 0;
    //task b
    //double rho_max = 1;
    //task c
    double rho_max = 4.5;
    double h = (rho_max-rho0)/n;
    double a = -1/(h*h);
    double d = 2/(h*h);

    //task b, potential = false
    //A = jacobi_method.create_matrix(A,n, d, a, rho_max, false, 0);
    //task c, potential = true

    A = jacobi_method.create_matrix(A,n, d, a, rho_max, true, 0.0);
    //task d, potential = true, omega != 0.0
    //A = jacobi_method.create_matrix(A,n, d, a, rho_max, true, 0.01);

    vec eig = eig_sym(A);
    while ( maxnondiag > tolerance && iterations <= maxiter)
    {
       int p, q;
       maxnondiag = jacobi_method.offdiag(A, p, q, n);
       jacobi_method.Jacobi_rotate(A, R, p, q, n);
       iterations++;
       /*if (iterations ==10000 or iterations == 15000 or iterations == 200000 or iterations == 300000 or iterations == 350000 ){
           cout << iterations << endl;
       }
        */

    }


    cout << "iterations: "<<iterations << endl;
    vec Eigenvectors = diagvec(A);
    vec eigenvec_sort = sort(Eigenvectors);
    //vec eig_sort = sort(eig);
    cout << eigenvec_sort(0) << " " << eigenvec_sort(1) <<" "  <<eigenvec_sort(2) << " "<<eigenvec_sort(3) << endl;
    //cout << eig_sort(0) << " " << eig_sort(1) << " " << eig_sort(2) <<eig_sort(3) << " " << eig_sort(4) << endl;

    return 0;
}



