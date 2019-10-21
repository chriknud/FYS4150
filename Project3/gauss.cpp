#include "gauss.h"
#include <mpi.h>

gauss::gauss()
{

}

double gauss::gammln( double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}
void gauss::gauss_laguerre(double *x, double *w, int n, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=n;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
}

void gauss::gauleg(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
       ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
     p2 =0.0;

       /*
       ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

     for(j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
     }

       /*
       ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */

     pp = n * (z * p1 - p2)/(z * z - 1.0);
     z1 = z;
     z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /*
      ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()

double gauss::func_cartesian(double x1, double y1, double z1, double x2, double y2, double z2){
    if  ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) > ZERO)
        return exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2)))
                  / sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    else
        return 0;
}
double gauss::func_polar(double r1, double theta1, double phi1, double r2, double theta2, double phi2){
    double cosb = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
    double f = exp(-3*(r1+r2))*r1*r1*r2*r2*sin(theta1)*sin(theta2)/sqrt(r1*r1+r2*r2-2*r1*r2*cosb);
    if(r1*r1+r2*r2-2*r1*r2*cosb > ZERO)
        return f;
    else
        return 0;
}

void gauss::monte_carlo(int n, double a, double b){
     double x1, y1, z1, x2, y2, z2;
     double MCint, MCintsqr2, fx, Variance;
     MCint = MCintsqr2=0.;
     double invers_period = 1./RAND_MAX;
     srand(time(NULL));
     for ( int i = 1;  i <= n; i++){
           x1 = a + (b-a)*double(rand())*invers_period;
           y1 = a + (b-a)*double(rand())*invers_period;
           z1 = a + (b-a)*double(rand())*invers_period;
           x2 = a + (b-a)*double(rand())*invers_period;
           y2 = a + (b-a)*double(rand())*invers_period;
           z2 = a + (b-a)*double(rand())*invers_period;
           fx = func_cartesian(x1, y1, z1, x2, y2, z2);
           MCint += fx;
           MCintsqr2 += fx*fx;
     }
     MCint = MCint/((double) n );
     MCintsqr2 = MCintsqr2/((double) n );
     double variance=(MCintsqr2-MCint*MCint)*pow(b-a, 6);
     MCint = MCint*pow(b-a, 6);
     double exact ;
     double pi=3.141592653589793238462643383279502884197;
     exact = 5*pi*pi/(16*16);
     cout << "Monte Carlo"<< endl;
     cout << "Number of integration points " << n << endl;
     cout << "Variance= " << variance << " Integral = " << MCint << " Exact= " << exact << endl;
     cout << "Rel. error = " << (abs(exact-MCint))/exact << endl;
    }

void gauss::monte_carlo_improved(int n){
    double r1, theta1, phi1, r2, theta2, phi2;
    double MCint, MCintsqr2, fx, Variance;
    double pi=3.141592653589793238462643383279502884197;
    MCint = MCintsqr2=0.;
    double invers_period = 1./RAND_MAX;
    srand(time(NULL));
    for ( int i = 1;  i <= n; i++){
          r1 = -log(1-double(rand())*invers_period);
          theta1 = pi*double(rand())*invers_period;
          phi1 = 2*pi*double(rand())*invers_period;
          r2 = -log(1-double(rand())*invers_period);
          theta2 = pi*double(rand())*invers_period;
          phi2 = 2*pi*double(rand())*invers_period;
          fx = func_polar(r1, theta1, phi1, r2, theta2, phi2);
          MCint += fx;
          MCintsqr2 += fx*fx;
    }
    MCint = MCint/((double) n );
    MCintsqr2 = MCintsqr2/((double) n );
    double variance=(MCintsqr2-MCint*MCint)*4*pow(pi, 4);
    MCint = MCint*4*pow(pi, 4);
    double exact ;
    exact = 5*pi*pi/(16*16);
    //tester_func(MCint, exact);
    cout << "Monte Carlo Improved"<< endl;
    cout << "Number of integration points " << n << endl;
    cout << "Variance= " << variance << " Integral = " << MCint << " Exact= " << exact << endl;
    cout << "Rel. error = " << (abs(exact-MCint))/exact << endl;
}

void gauss::monte_carlo_improved_MPI(int n){
    int local_n, numprocs, my_rank;
    double time_start, time_end, total_time;
    double total_sum;
    MPI_Init (NULL,NULL);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    time_start = MPI_Wtime();
    local_n = n/numprocs;
    double r1, theta1, phi1, r2, theta2, phi2;
    double MCint, MCintsqr2, fx, Variance;
    double pi=3.141592653589793238462643383279502884197;
    double exact = 5*pi*pi/(16*16);
    MCint = MCintsqr2 = 0.;
    double invers_period = 1./RAND_MAX;
    srand(time(NULL));
    for ( int i = 1;  i <= local_n; i++){
          r1 = -log(1-double(rand())*invers_period);
          theta1 = pi*double(rand())*invers_period;
          phi1 = 2*pi*double(rand())*invers_period;
          r2 = -log(1-double(rand())*invers_period);
          theta2 = pi*double(rand())*invers_period;
          phi2 = 2*pi*double(rand())*invers_period;
          fx = func_polar(r1, theta1, phi1, r2, theta2, phi2);
          MCint += fx;
          MCintsqr2 += fx*fx;
    }
    MCint = MCint/((double) n );
    MCintsqr2 = MCintsqr2/((double) n );
    double variance=(MCintsqr2-MCint*MCint)*4*pow(pi, 4);
    MCint = MCint*4*pow(pi, 4);
    MPI_Reduce(&MCint, &total_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    time_end = MPI_Wtime();
    cout << "hei" << my_rank << endl;
    if (my_rank==0){
        cout << "Monte Carlo Improved MPI"<< endl;
        cout << "Number of integration points " << n << endl;
        cout << "Variance = " << variance << " Integral = " << total_sum << " Exact= " << exact << endl;
        cout << "Time: " << time_end-time_start<< endl;
        cout << "Rel. error = " << (abs(exact-MCint))/exact << endl;
    }
    MPI_Finalize();
}


void gauss::tester_func(double final_MCint, double exact){
    double tolerance = 1e-3;
    if (abs(exact - final_MCint)<tolerance){
        cout << "Test for checking the exact value of the integral is passed." << endl;
    }
    else{
        cout << "Test not passed, the calculated values and the exact values are not the same." << endl;
        cout << "Integral    = " << final_MCint << endl;
        cout << "Exact       = " << exact << endl;
    }
}
