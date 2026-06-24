#include "jacobi.h"
#include <stdexcept>

namespace pp {

// A <- A * J
// This rotates columns p and q.
void Jacobi::timesJ(matrix& A, int p, int q, double theta)
{
    const double c = std::cos(theta);
    const double s = std::sin(theta);

    for(int i = 0; i < A.size1(); i++) {
        const double aip = A(i,p);
        const double aiq = A(i,q);

        A(i,p) = c*aip - s*aiq;
        A(i,q) = s*aip + c*aiq;
    }
}

// A <- J * A
// This rotates rows p and q.
void Jacobi::Jtimes(matrix& A, int p, int q, double theta)
{
    const double c = std::cos(theta);
    const double s = std::sin(theta);

    for(int j = 0; j < A.size2(); j++) {
        const double apj = A(p,j);
        const double aqj = A(q,j);

        A(p,j) =  c*apj + s*aqj;
        A(q,j) = -s*apj + c*aqj;
    }
}

// Cyclic Jacobi diagonalization.
// Input: symmetric matrix A.
// Output: eigenvalues in w, eigenvectors in columns of V.
int Jacobi::cyclic(matrix& A, vector& w, matrix& V, double tolerance)
{
    const int n = A.size1();

    if(A.size2() != n) {
        throw std::invalid_argument("Jacobi::cyclic: matrix must be square");
    }

    w.resize(n);
    V.resize(n,n);
    V.setid();

    bool changed;
    int sweeps = 0;

    do {
        changed = false;
        sweeps++;

        for(int p = 0; p < n-1; p++) {
            for(int q = p+1; q < n; q++) {
                const double apq = A(p,q);
                const double app = A(p,p);
                const double aqq = A(q,q);

                const double theta = 0.5 * std::atan2(2*apq, aqq-app);

                const double c = std::cos(theta);
                const double s = std::sin(theta);

                const double new_app = c*c*app - 2*s*c*apq + s*s*aqq;
                const double new_aqq = s*s*app + 2*s*c*apq + c*c*aqq;

                if(!approx(new_app, app, tolerance, tolerance) ||
                   !approx(new_aqq, aqq, tolerance, tolerance)) {

                    changed = true;

                    // A <- J^T A J
                    timesJ(A, p, q, theta);
                    Jtimes(A, p, q, -theta);

                    // V <- V J
                    timesJ(V, p, q, theta);
                }
            }
        }

    } while(changed);

    for(int i = 0; i < n; i++) {
        w[i] = A(i,i);
    }

    return sweeps;
}

// Constructor.
// Takes A by value, so the original matrix is not destroyed.
Jacobi::Jacobi(matrix A, double tolerance)
{
    cyclic(A, w, V, tolerance);
}

} // namespace pp