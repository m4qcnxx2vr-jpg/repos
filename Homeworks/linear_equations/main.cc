#include "matrix.h"
#include "qr.h"
#include <string>
#include <random>
#include <iostream>
#include <cmath>
#include <cassert>

// ------------------------------------------------------------
// Helper: identity matrix
// ------------------------------------------------------------
pp::matrix eye(int n){
    pp::matrix M(n,n,0.0);
    for(int i=0;i<n;i++) M(i,i)=1;
    return M;
}

// ------------------------------------------------------------
// Helper: random vector
// ------------------------------------------------------------
pp::vector random_vector(int n, unsigned seed = 2){
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0,1.0);

    pp::vector v(n);
    for(int i=0;i<n;i++) v[i] = dist(rng);
    return v;
}

// ------------------------------------------------------------
// Helper: random matrix (square or rectangular)
// ------------------------------------------------------------
pp::matrix random_matrix(int n, int m, unsigned seed = 1){
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0,1.0);

    pp::matrix A(n,m,0.0);
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            A(i,j) = dist(rng);

    return A;
}

// ------------------------------------------------------------
// Helper: random tall matrix A (n > m)
// ------------------------------------------------------------
pp::matrix random_tall(int n, int m, unsigned seed = 1){
    assert(n > m);
    return random_matrix(n, m, seed);
}

// ------------------------------------------------------------
// Helper: check if R is upper triangular
// ------------------------------------------------------------
bool is_upper_triangular(const pp::matrix& R, double tol = 1e-10){
    if(R.size1() != R.size2()) return false;
    int m = R.size1();
    for(int i=0;i<m;i++)
        for(int j=0;j<i;j++) // below diagonal
            if(std::fabs(R(i,j)) > tol) return false;
    return true;
}

// ------------------------------------------------------------
// Helper: matrix approx
// ------------------------------------------------------------
bool approx_mat(const pp::matrix& A, const pp::matrix& B, double tol = 1e-6){
    if(A.size1() != B.size1() || A.size2() != B.size2()) return false;
    for(int i=0;i<A.size1();i++)
        for(int j=0;j<A.size2();j++)
            if(std::fabs(A(i,j) - B(i,j)) > tol) return false;
    return true;
}

int main(int argc, char** argv){


    // ------------------------------------------------------------
// Benchmark mode: ./main -size N
// Only do QR factorization on a random NxN matrix and exit.
// ------------------------------------------------------------
if(argc == 3 && std::string(argv[1]) == "-size") {
    int N = std::stoi(argv[2]);

    pp::matrix A = random_matrix(N, N, 123); // random NxN
    pp::qr decomp(A);                        // factorize

    // optional: print something tiny so compiler doesn't optimize away (it won't)
    // std::cout << decomp.R(0,0) << "\n";

        return 0;
    }
    // ========================================================
    // VECTOR TESTS
    // ========================================================
    std::cout << "====================================\n";
    std::cout << "VECTOR TESTS\n";
    std::cout << "====================================\n";
    std::cout << "Testing vector construction, copying, and arithmetic.\n\n";

    int n = 3;

    pp::vector v(n);
    std::cout << "Initial vector (should be zeros):\n";
    v.print("v=");

    std::cout << "\nAssign v[i] = i+1:\n";
    for(int i=0;i<v.size();i++) v[i]=i+1;
    v.print("v=");

    std::cout << "\nCopy v into u:\n";
    pp::vector u = v;
    u.print("u=");

    std::cout << "\nModify v[0] = 999. u should stay unchanged:\n";
    v[0] = 999;
    v.print("v=");
    u.print("u=");

    std::cout << "\nCreate w as a copy of v, then do w += v:\n";
    pp::vector w(v);
    w += v;
    w.print("w+=v");

    std::cout << "\nTest scalar division (w/10):\n";
    (w/10).print("w/10");

    std::cout << "\nTest vector addition (u+u):\n";
    (u+u).print("u+u");


    // ========================================================
    // MATRIX TESTS
    // ========================================================
    std::cout << "\n====================================\n";
    std::cout << "MATRIX TESTS\n";
    std::cout << "====================================\n";
    std::cout << "Testing matrix creation, arithmetic, and multiplication.\n\n";

    pp::matrix M(n,n);

    std::cout << "Initial matrix M (should be zeros):\n";
    M.print("M=");

    std::cout << "\nFill M(i,j) = i + j:\n";
    for(int i=0;i<M.size1();i++)
        for(int j=0;j<M.size2();j++)
            M(i,j) = i + j;

    M.print("M=");

    std::cout << "\nTest matrix addition (M+M):\n";
    (M+M).print("M+M");

    std::cout << "\nTest matrix subtraction (M-M):\n";
    (M-M).print("M-M");

    std::cout << "\nTest chained multiplication: M*(M+M)*M\n";
    (M*(M+M)*M).print("M*(M+M)*M");

    std::cout << "\nTest matrix-vector multiplication: M*u\n";
    (M*u).print("M*u");

    std::cout << "\nTest scalar division of matrix: M/10\n";
    (M/10).print("M/10");

    std::cout << "\nIdentity matrix I:\n";
    pp::matrix I3 = eye(n);
    I3.print("I=");


    // ========================================================
    // QR DECOMPOSITION TESTS (tall matrix)
    // ========================================================
    {
        std::cout << "\n====================================\n";
        std::cout << "QR DECOMPOSITION TESTS\n";
        std::cout << "====================================\n";
        std::cout << "Checks required by the assignment:\n";
        std::cout << "  1) R is upper triangular\n";
        std::cout << "  2) Q^T Q = I\n";
        std::cout << "  3) Q R = A\n\n";

        int rows = 8;
        int cols = 5;

        std::cout << "Generating random tall matrix A (" << rows << " x " << cols << ")...\n";
        pp::matrix A = random_tall(rows, cols, 42);
        A.print("\nMatrix A:");

        std::cout << "\nFactorizing A = Q*R using stabilized Gram-Schmidt...\n";
        pp::qr decomp(A);

        pp::matrix Q = decomp.Q;
        pp::matrix R = decomp.R;

        Q.print("\nMatrix Q (orthonormal columns):");
        R.print("\nMatrix R (upper triangular):");

        std::cout << "\n(1) Checking if R is upper triangular...\n";
        std::cout << "Result: " << (is_upper_triangular(R) ? "OK" : "FAIL") << "\n";

        pp::matrix QtQ = Q.transpose() * Q;
        QtQ.print("\nMatrix Q^T Q:");

        pp::matrix Icols = eye(cols);
        Icols.print("\nIdentity matrix I:");

        std::cout << "\n(2) Checking orthogonality: Q^T Q = I...\n";
        std::cout << "Result: " << (approx_mat(QtQ, Icols, 1e-6) ? "OK" : "FAIL") << "\n";

        pp::matrix QR = Q * R;
        QR.print("\nMatrix Q*R:");

        std::cout << "\n(3) Checking reconstruction: Q*R = A...\n";
        std::cout << "Result: " << (approx_mat(QR, A, 1e-6) ? "OK" : "FAIL") << "\n";

        std::cout << "\n====================================\n";
        std::cout << "END OF QR TESTS\n";
        std::cout << "====================================\n";
    }


    // ========================================================
    // QR SOLVE TEST (square matrix)
    // ========================================================
    {
        std::cout << "\n====================================\n";
        std::cout << "QR SOLVE TEST\n";
        std::cout << "====================================\n";

        int nsolve = 6;

        std::cout << "Generating random square matrix A (" << nsolve << " x " << nsolve << ")...\n";
        pp::matrix A2 = random_matrix(nsolve, nsolve, 5);
        A2.print("\nMatrix A:");

        std::cout << "\nGenerating random vector b...\n";
        pp::vector b = random_vector(nsolve, 7);
        b.print("Vector b:");

        std::cout << "\nFactorizing A = Q*R...\n";
        pp::qr decomp2(A2);

        std::cout << "\nSolving Ax = b using QR...\n";
        pp::vector x = decomp2.solve(b);
        x.print("Computed solution x:");

        std::cout << "\nChecking Ax...\n";
        pp::vector Ax = A2 * x;
        Ax.print("A*x:");
        b.print("Original b:");

        std::cout << "\nChecking if A*x ≈ b...\n";
        bool ok = true;
        for(int i=0;i<nsolve;i++)
            if(!pp::approx(Ax[i], b[i], 1e-6, 1e-6))
                ok = false;

        std::cout << "Result: " << (ok ? "OK" : "FAIL") << "\n";

        std::cout << "\n====================================\n";
        std::cout << "END OF SOLVE TEST\n";
        std::cout << "====================================\n";
    }


    // ========================================================
    // DETERMINANT TEST
    // ========================================================
    {
        std::cout << "\n====================================\n";
        std::cout << "DETERMINANT TEST\n";
        std::cout << "====================================\n";

        int ndet = 4;
        pp::matrix A3 = random_matrix(ndet, ndet, 10);
        A3.print("Matrix A:");

        pp::qr decomp3(A3);
        double d = decomp3.det();

        std::cout << "\nDeterminant det(A) = " << d << "\n";
        std::cout << "\n====================================\n";
        std::cout << "END OF DETERMINANT TEST\n";
        std::cout << "====================================\n";
    }


    // ========================================================
    // INVERSE TEST (via QR)
    // ========================================================
    {
        std::cout << "\n====================================\n";
        std::cout << "INVERSE TEST (via QR)\n";
        std::cout << "====================================\n";
        std::cout << "We check that A * A^{-1} = I.\n\n";

        int ninv = 5;

        std::cout << "Generating random square matrix A (" << ninv << " x " << ninv << ")...\n";
        pp::matrix Ainv_test = random_matrix(ninv, ninv, 123);
        Ainv_test.print("\nMatrix A:");

        std::cout << "\nFactorizing A = Q*R ...\n";
        pp::qr decomp_inv(Ainv_test);

        std::cout << "\nComputing inverse B = A^{-1} ...\n";
        pp::matrix B = decomp_inv.inverse();
        B.print("\nMatrix B = A^{-1}:");

        std::cout << "\nComputing A*B ... (should be identity)\n";
        pp::matrix AB = Ainv_test * B;
        AB.print("\nMatrix A*B:");

        pp::matrix Iinv = eye(ninv);
        Iinv.print("\nIdentity I:");

        std::cout << "\nChecking if A*B ≈ I ...\n";
        std::cout << "Result: " << (approx_mat(AB, Iinv, 1e-6) ? "OK" : "FAIL") << "\n";

        std::cout << "\n====================================\n";
        std::cout << "END OF INVERSE TEST\n";
        std::cout << "====================================\n";
    }

    return 0;
}