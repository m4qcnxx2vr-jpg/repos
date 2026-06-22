#include "matrix.h"
#include "svd.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

double matrix_norm(const pp::matrix& A) {
    double sum = 0.0;

    for (int i = 0; i < A.size1(); ++i) {
        for (int j = 0; j < A.size2(); ++j) {
            sum += A(i, j) * A(i, j);
        }
    }

    return std::sqrt(sum);
}

double offdiag_norm(const pp::matrix& A) {
    double sum = 0.0;

    for (int i = 0; i < A.size1(); ++i) {
        for (int j = 0; j < A.size2(); ++j) {
            if (i != j) {
                sum += A(i, j) * A(i, j);
            }
        }
    }

    return std::sqrt(sum);
}

bool approx_matrix(const pp::matrix& A, const pp::matrix& B,
                   double acc = 1e-6, double eps = 1e-6) {
    if (A.size1() != B.size1() || A.size2() != B.size2()) {
        return false;
    }

    for (int i = 0; i < A.size1(); ++i) {
        for (int j = 0; j < A.size2(); ++j) {
            double diff = std::abs(A(i, j) - B(i, j));
            double scale = std::max(std::abs(A(i, j)), std::abs(B(i, j)));

            if (diff > acc + eps * scale) {
                return false;
            }
        }
    }

    return true;
}

pp::matrix make_test_matrix(int n) {
    pp::matrix A(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A(i, j) = 1.0 / (1.0 + i + j);

            if (i == j) {
                A(i, j) += n;
            }
        }
    }

    return A;
}

void test_svd(const pp::matrix& A, const std::string& name) {
    std::cout << "\n----------------------------------------\n";
    std::cout << name << "\n";
    std::cout << "----------------------------------------\n";

    pp::SVD svd(A);

    pp::matrix reconstructed = svd.U * svd.D * svd.V.T();

    pp::matrix IU(svd.U.size2(), svd.U.size2());
    pp::matrix IV(svd.V.size2(), svd.V.size2());

    IU.setid();
    IV.setid();

    pp::matrix UtU = svd.U.T() * svd.U;
    pp::matrix VtV = svd.V.T() * svd.V;

    std::cout << "A =\n";
    A.print();

    std::cout << "\nU =\n";
    svd.U.print();

    std::cout << "\nD =\n";
    svd.D.print();

    std::cout << "\nV =\n";
    svd.V.print();

    std::cout << "\nSingular values =\n";
    svd.s.print();

    std::cout << "\nU D V^T =\n";
    reconstructed.print();

    std::cout << "\nU^T U =\n";
    UtU.print();

    std::cout << "\nV^T V =\n";
    VtV.print();

    std::cout << "\nTests:\n";
    std::cout << "U D V^T approx A : "
              << approx_matrix(reconstructed, A) << "\n";

    std::cout << "U^T U approx I   : "
              << approx_matrix(UtU, IU) << "\n";

    std::cout << "V^T V approx I   : "
              << approx_matrix(VtV, IV) << "\n";

    std::cout << "D offdiag norm   : "
              << offdiag_norm(svd.D) << "\n";

    std::cout << "reconstruction error norm : "
              << matrix_norm(reconstructed - A) << "\n";
}

void benchmark_svd(int n) {
    pp::matrix A = make_test_matrix(n);
    pp::SVD svd(A);
}

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "-size") {
            int n = std::stoi(argv[2]);
            benchmark_svd(n);
            return 0;
        }

        if (argc != 1) {
            std::cerr << "Usage: " << argv[0] << " [-size N]\n";
            return 1;
        }

        pp::matrix A(3, 3);

        A(0, 0) = 1.0;  A(0, 1) = 2.0;  A(0, 2) = 3.0;
        A(1, 0) = 4.0;  A(1, 1) = 5.0;  A(1, 2) = 6.0;
        A(2, 0) = 7.0;  A(2, 1) = 8.0;  A(2, 2) = 10.0;

        test_svd(A, "Test 1: general 3x3 matrix");

        pp::matrix B(3, 3);

        B(0, 0) = 3.0;  B(0, 1) = 0.0;  B(0, 2) = 0.0;
        B(1, 0) = 0.0;  B(1, 1) = 2.0;  B(1, 2) = 0.0;
        B(2, 0) = 0.0;  B(2, 1) = 0.0;  B(2, 2) = 1.0;

        test_svd(B, "Test 2: diagonal matrix");

        pp::matrix C(2, 2);

        C(0, 0) = 1.0;  C(0, 1) = 2.0;
        C(1, 0) = 3.0;  C(1, 1) = 4.0;

        test_svd(C, "Test 3: general 2x2 matrix");

        pp::matrix E(3, 3);

        E(0, 0) = -1.0;  E(0, 1) =  2.0;  E(0, 2) = -3.0;
        E(1, 0) =  4.0;  E(1, 1) = -5.0;  E(1, 2) =  6.0;
        E(2, 0) = -7.0;  E(2, 1) =  8.0;  E(2, 2) = -9.0;

        test_svd(E, "Test 4: matrix with negative entries");

    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }

    return 0;
}