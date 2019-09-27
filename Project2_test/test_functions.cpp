#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "jacobi.h"
jacobi jacobi_method;

TEST_CASE("Testing max A(i,j)"){
    int n = 4;
    mat A(n,n,fill::zeros);

    double d = 2.0; //(h*h);
    double a = -1.0; //(h*h);
    int p, q;
    A = jacobi_method.create_matrix(A, n, d, a);
    REQUIRE(A(2,2)==2);
    REQUIRE(A(1,2)==-1);
    A(1,3) = 5.6;
    double maxnondiag = jacobi_method.offdiag(A, p, q, n);
    REQUIRE(maxnondiag==Approx(5.6));
}

