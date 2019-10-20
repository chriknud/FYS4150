#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "gauss.h"
#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-10
using namespace std;

void loopyloop_cart(int N, double a, double b);
void loopyloop_polar(int N);
int main()
{
    gauss gauss;
     int n;
     double a, b;
     n = 100000;
     a = -3;
     b = 3;
/*
//Loop over Legendre
     for (int i=1; i<=7;i++){
        int N = i *5;
        auto t1a = std::chrono::high_resolution_clock::now();
        loopyloop_cart(N, a, b);
        auto t2a = std::chrono::high_resolution_clock::now();
        auto duration_legendre = std::chrono::duration_cast<std::chrono::microseconds>( t2a - t1a ).count();
        cout << "Time used with Gauss-Legendre quadrature: " << duration_legendre*1e-6 << " sec." << endl;

     }

//Loop over Laguerre
     for (int i=1; i<=7;i++){
        int N = i *5;
        auto t1b = std::chrono::high_resolution_clock::now();
        loopyloop_polar(N);
        auto t2b = std::chrono::high_resolution_clock::now();
        auto duration_laguerre = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();
        cout << "Time used with Gauss-Laguerre quadrature: " << duration_laguerre*1e-6 << " sec." << endl;

     }
//Loop over Monte Carlo
     for (int i=5; i<=8;i++){
        double N = pow(10,i);
        int numb = N;
        cout << "N= " << N << endl;
        auto t1c = std::chrono::high_resolution_clock::now();
        gauss.monte_carlo(numb, a, b);
        auto t2c = std::chrono::high_resolution_clock::now();
        auto duration_monte_carlo = std::chrono::duration_cast<std::chrono::microseconds>( t2c - t1c ).count();
        cout << "Time used with Monte Carlo: " << duration_monte_carlo*1e-6 << " sec." << endl;
     }
//Loop over Monte Carlo improved
     for (int i=5; i<=8;i++){
        double N = pow(10,i);
        int numb = N;
        cout << "N= " << N << endl;
        auto t1d = std::chrono::high_resolution_clock::now();
        gauss.monte_carlo_improved(numb, a, b);
        auto t2d = std::chrono::high_resolution_clock::now();
        auto duration_monte_carlo_improved = std::chrono::duration_cast<std::chrono::microseconds>( t2d - t1d ).count();
        cout << "Time used with improved Monte Carlo: " << duration_monte_carlo_improved*1e-6 << " sec." << endl;
     }
     */

     int N1 = 25;
     int N2 = 15;

//Gauss Legendre quadrature
     auto t1a = std::chrono::high_resolution_clock::now();
     loopyloop_cart(N1, a, b);
     auto t2a = std::chrono::high_resolution_clock::now();
     auto duration_legendre = std::chrono::duration_cast<std::chrono::microseconds>( t2a - t1a ).count();
     cout << "Time used with Gauss-Legendre quadrature: " << duration_legendre*1e-6 << " sec." << endl;

//Gauss Laguerre quadrature
     auto t1b = std::chrono::high_resolution_clock::now();
     loopyloop_polar(N2);
     auto t2b = std::chrono::high_resolution_clock::now();
     auto duration_laguerre = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();
     cout << "Time used with Gauss-Laguerre quadrature: " << duration_laguerre*1e-6 << " sec." << endl;

//Monte Carlo
     auto t1c = std::chrono::high_resolution_clock::now();
     gauss.monte_carlo(n, a, b);
     auto t2c = std::chrono::high_resolution_clock::now();
     auto duration_monte_carlo = std::chrono::duration_cast<std::chrono::microseconds>( t2c - t1c ).count();
     cout << "Time used with Monte Carlo: " << duration_monte_carlo*1e-6 << " sec." << endl;

//Monte Carlo improved
     auto t1d = std::chrono::high_resolution_clock::now();
     gauss.monte_carlo_improved(n, a, b);
     auto t2d = std::chrono::high_resolution_clock::now();
     auto duration_monte_carlo_improved = std::chrono::duration_cast<std::chrono::microseconds>( t2d - t1d ).count();
     cout << "Time used with improved Monte Carlo: " << duration_monte_carlo_improved*1e-6 << " sec." << endl;

//Monte Carlo improved with MPI
     gauss.monte_carlo_improved_MPI(n, a, b);

}

void loopyloop_cart(int N, double a, double b){
    double exact ;
    double pi=3.141592653589793238462643383279502884197;
    exact = 5*pi*pi/(16*16);
    gauss gauss;
    double *x = new double [N];
    double *w = new double [N];
    gauss.gauleg(a,b,x,w, N);
    double int_gauss = 0.;
    for (int i=0;i<N;i++){
        for (int j = 0;j<N;j++){
        for (int k = 0;k<N;k++){
        for (int l = 0;l<N;l++){
        for (int m = 0;m<N;m++){
        for (int n = 0;n<N;n++){
       int_gauss+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]
      *gauss.func_cartesian(x[i],x[j],x[k],x[l],x[m],x[n]);
           }}}}}
       }
        cout << "Gauss Legendre Quadrature"<< endl;
        cout << "Number of integration points " << N << endl;
        cout << "Exact value " << exact << endl;
        cout << "Numerical value " << int_gauss << endl;

}

void loopyloop_polar(int N){
    double exact ;
    double pi=3.141592653589793238462643383279502884197;
    exact = 5*pi*pi/(16*16);
    gauss gauss;
    double *theta = new double [N];
    double *phi = new double [N];
    double *r = new double [N];
    double *w_t = new double[N];
    double *w_p = new double [N];
    double *w_r = new double [N];
    gauss.gauleg(0,pi,theta, w_t, N);
    gauss.gauleg(0,2*pi,phi, w_p, N);
    gauss.gauss_laguerre(r, w_r, N, 0);
    double int_gauss = 0.;
    for (int i=0;i<N;i++){
        for (int j = 0;j<N;j++){
        for (int k = 0;k<N;k++){
        for (int l = 0;l<N;l++){
        for (int m = 0;m<N;m++){
        for (int n = 0;n<N;n++){
       int_gauss+=w_r[i+1]*w_t[j]*w_p[k]*w_r[l+1]*w_t[m]*w_p[n]
      *gauss.func_polar(r[i+1],theta[j],phi[k],r[l+1],theta[m],phi[n]);
           }}}}}
   }
    cout << "Gauss Laguerre Quadrature"<< endl;
    cout << "Number of integration points " << N << endl;
    cout << "Exact value " << exact << endl;
    cout << "Numerical value " << int_gauss << endl;
}

