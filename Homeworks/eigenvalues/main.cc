#include "matrix.h"
#include "jacobi.h"
#include <iostream>
#include <cmath>
#include <algorithm>

#include "matrix.h"
#include "jacobi.h"
#include <random>
#include <iostream>

static pp::matrix diag_from_w(const pp::vector& w){
    int n = w.size();
    pp::matrix D(n,n,0.0);
    for(int i=0;i<n;i++) D(i,i) = w[i];
    return D;
}


int main(){

    int n = 5;

    std::mt19937 rng(1);
    std::uniform_real_distribution<double> U(-1,1);

    // random symmetric matrix
    pp::matrix A(n,n);
    for(int i=0;i<n;i++)
        for(int j=i;j<n;j++){
            double x = U(rng);
            A(i,j) = x;
            A(j,i) = x;
        }

    A.print("Random symmetric matrix A:");

    pp::Jacobi evd(A);

    pp::matrix D = diag_from_w(evd.w);

    pp::matrix VTAV = evd.V.T() * A * evd.V;
    pp::matrix VDVT = evd.V * D * evd.V.T();
    pp::matrix VTV  = evd.V.T() * evd.V;
    pp::matrix VVT  = evd.V * evd.V.T();

    VTAV.print("V^T A V (should equal D):");
    D.print("D:");

    VDVT.print("V D V^T (should equal A):");
    A.print("A:");

    VTV.print("V^T V (should be identity):");
    VVT.print("V V^T (should be identity):");

    return 0;
}