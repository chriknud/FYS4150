#ifndef JACOBI_H
#define JACOBI_H
#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include "armadillo"
using namespace std;
using namespace arma;

class jacobi
{
public:
    jacobi();
    double offdiag(mat& A, int& p, int& q, int n);
    void Jacobi_rotate ( mat &A, mat &R, int k, int l, int n );
    mat create_matrix(mat M,int n, double d, double a, double rho_max, bool potential, double omega);
};

#endif // JACOBI_H
