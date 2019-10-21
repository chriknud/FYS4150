#ifndef GAUSS_H
#define GAUSS_H
#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include "armadillo"
#include <mpi.h>
using namespace std;
using namespace arma;
#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-10

class gauss
{
public:
    gauss();
    void gauss_laguerre(double *x, double *w, int n, double alf);
    void gauleg(double x1, double x2, double x[], double w[], int n);
    double func_cartesian(double x1, double y1, double z1, double x2, double y2, double z2);
    double gammln( double xx);
    double func_polar(double r1, double theta1, double phi1, double r2, double theta2, double phi2);
    void monte_carlo(int n, double a, double b);
    void monte_carlo_improved(int n, double a, double b);
    void monte_carlo_improved_MPI(int n, double a, double b);
    void tester_func(double final_MCint, double exact);

};

#endif // GAUSS_H
