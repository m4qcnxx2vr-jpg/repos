#include "svd.h"
#include <stdexcept>
#include <cmath>

namespace pp {

SVD::SVD(matrix A, double tolerance) {
    if (A.size1() != A.size2()) {
        throw std::invalid_argument("SVD: matrix must be square");
    }

    int n = A.size1();

    U = matrix(n, n);
    V = matrix(n, n);
    D = A;
    s = vector(n);

    U.setid();
    V.setid();

    two_sided(D, U, V, s, tolerance);
}

void SVD::timesJ(matrix& A, int p, int q, double theta) {
    double c = std::cos(theta);
    double s = std::sin(theta);

    for (int i = 0; i < A.size1(); ++i) {
        double Aip = A(i, p);
        double Aiq = A(i, q);

        A(i, p) = c * Aip - s * Aiq;
        A(i, q) = s * Aip + c * Aiq;
    }
}

void SVD::Jtimes(matrix& A, int p, int q, double theta) {
    double c = std::cos(theta);
    double s = std::sin(theta);

    for (int j = 0; j < A.size2(); ++j) {
        double Apj = A(p, j);
        double Aqj = A(q, j);

        A(p, j) = c * Apj - s * Aqj;
        A(q, j) = s * Apj + c * Aqj;
    }
}

int SVD::two_sided(matrix& A, matrix& U, matrix& V, vector& sing,
                   double tolerance) {
    if (A.size1() != A.size2()) {
        throw std::invalid_argument("SVD::two_sided: matrix must be square");
    }

    int n = A.size1();

    U = matrix(n, n);
    V = matrix(n, n);
    sing = vector(n);

    U.setid();
    V.setid();

    int sweeps = 0;
    int max_sweeps = 100 * n * n;

    bool changed = true;

    while (changed && sweeps < max_sweeps) {
        changed = false;
        ++sweeps;

        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {

                double app = A(p, p);
                double aqq = A(q, q);
                double apq = A(p, q);
                double aqp = A(q, p);

                double before = std::abs(apq) + std::abs(aqp);

                if (before <= tolerance) {
                    continue;
                }

                /*
                    First rotation:
                    A <- G^T A

                    The angle gamma is chosen such that after this
                    left rotation the two off-diagonal elements are equal:
                    A(p,q) = A(q,p).
                */
                double gamma = std::atan2(apq - aqp, aqq + app);

                Jtimes(A, p, q, gamma);   // A <- G^T A
                timesJ(U, p, q, gamma);   // U <- U G

                /*
                    Second rotation:
                    A <- J^T A J

                    Now the p,q block is symmetric, so we use the usual
                    Jacobi angle to remove the off-diagonal elements.
                */
                app = A(p, p);
                aqq = A(q, q);
                apq = A(p, q);

                double theta = 0.5 * std::atan2(2.0 * apq, aqq - app);

                Jtimes(A, p, q, theta);   // A <- J^T A
                timesJ(A, p, q, theta);   // A <- A J

                timesJ(U, p, q, theta);   // U <- U J
                timesJ(V, p, q, theta);   // V <- V J

                double after = std::abs(A(p, q)) + std::abs(A(q, p));

                if (after < before) {
                    changed = true;
                }
            }
        }
    }

    /*
        The diagonal entries of D should be non-negative.
        If D(i,i) is negative, multiply column i of U by -1.
    */
    for (int i = 0; i < n; ++i) {
        if (A(i, i) < 0.0) {
            A(i, i) = -A(i, i);

            for (int k = 0; k < n; ++k) {
                U(k, i) = -U(k, i);
            }
        }

        sing[i] = A(i, i);
    }

    /*
        Remove tiny numerical noise outside the diagonal.
    */
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j && std::abs(A(i, j)) < tolerance) {
                A(i, j) = 0.0;
            }
        }
    }

    return sweeps;
}

} 