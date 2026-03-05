#include "jacobi.h"
#include <algorithm>   // std::max
#include <stdexcept>   // std::invalid_argument

namespace pp {

// A <- A * J(p,q,theta)
// This rotates columns p and q of A.
void Jacobi::timesJ(matrix& A, int p, int q, double theta){
    const double c = std::cos(theta);
    const double s = std::sin(theta);

    for(int i=0;i<A.size1();i++){
        const double aip = A(i,p);
        const double aiq = A(i,q);
        A(i,p) = c*aip - s*aiq;
        A(i,q) = s*aip + c*aiq;
    }
}

// A <- J(p,q,theta)^T * A
// This rotates rows p and q of A.
void Jacobi::Jtimes(matrix& A, int p, int q, double theta){
    const double c = std::cos(theta);
    const double s = std::sin(theta);

    // Correct: A <- J^T * A, where J^T = [[c, -s],[s, c]]
    for(int j=0;j<A.size2();j++){
        const double apj = A(p,j);
        const double aqj = A(q,j);
        A(p,j) = c*apj - s*aqj;
        A(q,j) = s*apj + c*aqj;
    }
}

// Helper: maximum absolute off-diagonal element
static double max_offdiag(const matrix& A){
    const int n = A.size1();
    double m = 0.0;
    for(int i=0;i<n;i++)
        for(int j=i+1;j<n;j++)
            m = std::max(m, std::fabs(A(i,j)));
    return m;
}

// Cyclic Jacobi diagonalization (symmetric matrices only).
// On return: w holds diagonal (eigenvalues), V holds eigenvectors (columns).
int Jacobi::cyclic(matrix& A, vector& w, matrix& V, double tol){
    const int n = A.size1();
    if(A.size2()!=n) throw std::invalid_argument("Jacobi::cyclic: A must be square");

    w.resize(n);
    V.resize(n,n);
    V.setid();

    int sweeps = 0;

    while(true){
        sweeps++;

        // One sweep over all pairs (p,q)
        for(int p=0;p<n;p++){
            for(int q=p+1;q<n;q++){
                const double apq = A(p,q);
                if(std::fabs(apq) <= tol) continue;

                const double app = A(p,p);
                const double aqq = A(q,q);

                // theta chosen so the rotation kills A(p,q)
                // stable formula:
                const double theta = 0.5 * std::atan2(2*apq, (aqq - app));

                // Similarity transform: A <- J^T A J
                Jtimes(A, p, q, theta);
                timesJ(A, p, q, theta);

                // Accumulate eigenvectors: V <- V J
                timesJ(V, p, q, theta);
            }
        }

        // stop when off-diagonal elements are small
        if(max_offdiag(A) < tol) break;

        // optional safety stop
        if(sweeps > 50*n*n) break;
    }

    // eigenvalues = diagonal
    for(int i=0;i<n;i++) w[i] = A(i,i);

    return sweeps;
}

// Constructor: compute EVD of Ain
Jacobi::Jacobi(const matrix& Ain, double tol){
    matrix A = Ain;            // work copy (we diagonalize this)
    cyclic(A, w, V, tol);      // fills w and V
}

} // namespace pp