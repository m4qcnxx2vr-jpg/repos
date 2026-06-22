#include "jacobi.h"
#include "matrix.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using pp::matrix;
using pp::vector;

struct Options {
    double rmax = 10.0;
    double dr = 0.3;
    int nstates = 4;
};

Options read_options(int argc, char** argv)
{
    Options opt;

    for(int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if(arg == "-rmax" && i+1 < argc) {
            opt.rmax = std::atof(argv[++i]);
        }
        else if(arg == "-dr" && i+1 < argc) {
            opt.dr = std::atof(argv[++i]);
        }
        else if(arg == "-nstates" && i+1 < argc) {
            opt.nstates = std::atoi(argv[++i]);
        }
    }

    return opt;
}

matrix diagonal_matrix(const vector& w)
{
    matrix D(w.size(), w.size(), 0.0);

    for(int i = 0; i < w.size(); i++) {
        D(i,i) = w[i];
    }

    return D;
}

double max_abs(const matrix& A)
{
    double result = 0.0;

    for(int i = 0; i < A.size1(); i++) {
        for(int j = 0; j < A.size2(); j++) {
            result = std::max(result, std::fabs(A(i,j)));
        }
    }

    return result;
}

matrix identity(int n)
{
    matrix I(n,n,0.0);
    I.setid();
    return I;
}

matrix random_symmetric_matrix(int n)
{
    matrix A(n,n,0.0);

    for(int i = 0; i < n; i++) {
        for(int j = i; j < n; j++) {
            double x = std::sin(1.0 + 3.0*i + 7.0*j);

            A(i,j) = x;
            A(j,i) = x;
        }
    }

    return A;
}

matrix hydrogen_hamiltonian(double rmax, double dr)
{
    const int npoints = static_cast<int>(rmax/dr) - 1;

    matrix H(npoints, npoints, 0.0);

    const double factor = -0.5/(dr*dr);

    for(int i = 0; i < npoints; i++) {
        double r = dr*(i+1);

        H(i,i) = -2.0*factor - 1.0/r;

        if(i+1 < npoints) {
            H(i,i+1) = factor;
            H(i+1,i) = factor;
        }
    }

    return H;
}

void sort_eigenpairs(vector& w, matrix& V)
{
    const int n = w.size();

    for(int i = 0; i < n-1; i++) {
        int smallest = i;

        for(int j = i+1; j < n; j++) {
            if(w[j] < w[smallest]) {
                smallest = j;
            }
        }

        if(smallest != i) {
            std::swap(w[i], w[smallest]);

            for(int row = 0; row < n; row++) {
                std::swap(V(row,i), V(row,smallest));
            }
        }
    }
}

void part_A()
{
    std::cout << "A. Jacobi diagonalization with cyclic sweeps\n";

    matrix A = random_symmetric_matrix(5);
    matrix A_original = A;

    pp::Jacobi evd(A);

    vector w = evd.w;
    matrix V = evd.V;
    matrix D = diagonal_matrix(w);

    matrix VTAV = V.T()*A_original*V;
    matrix VDVT = V*D*V.T();
    matrix VTV = V.T()*V;
    matrix VVT = V*V.T();
    matrix I = identity(A.size1());

    std::cout << "max|V^T A V - D| = " << max_abs(VTAV - D) << '\n';
    std::cout << "max|V D V^T - A| = " << max_abs(VDVT - A_original) << '\n';
    std::cout << "max|V^T V - I|   = " << max_abs(VTV - I) << '\n';
    std::cout << "max|V V^T - I|   = " << max_abs(VVT - I) << '\n';

    std::cout << '\n';
}

void write_wavefunctions(double rmax, double dr, int nstates, const vector& w, const matrix& V)
{
    const int npoints = static_cast<int>(rmax/dr) - 1;
    const int states_to_write = std::min(nstates, w.size());

    const double norm_factor = 1.0/std::sqrt(dr);

    for(int k = 0; k < states_to_write; k++) {
        std::ofstream file("state_" + std::to_string(k) + ".dat");

        for(int i = 0; i < npoints; i++) {
            double r = dr*(i+1);
            double f = norm_factor*V(i,k);

            file << r << ' ' << f << '\n';
        }
    }
}

void part_B(double rmax, double dr, int nstates)
{
    std::cout << "B. Hydrogen atom, s-wave radial Schrödinger equation\n";
    std::cout << "Using rmax = " << rmax << ", dr = " << dr << '\n';

    matrix H = hydrogen_hamiltonian(rmax, dr);

    pp::Jacobi evd(H);

    vector w = evd.w;
    matrix V = evd.V;

    sort_eigenpairs(w, V);

    std::cout << "Lowest numerical energies compared with exact E_n = -1/(2 n^2):\n";
    std::cout << "state  numerical           exact               difference\n";

    const int states_to_print = std::min(nstates, w.size());

    for(int k = 0; k < states_to_print; k++) {
        double exact = -1.0/(2.0*(k+1)*(k+1));

        std::cout << std::setw(5) << k
                  << std::setw(20) << w[k]
                  << std::setw(20) << exact
                  << std::setw(20) << w[k] - exact << '\n';
    }

    write_wavefunctions(rmax, dr, nstates, w, V);

    std::cout << "Wavefunctions written to state_0.dat, state_1.dat, ...\n\n";

    std::cout << "Convergence with fixed rmax = " << rmax << " and changing dr:\n";
std::cout << "dr        ground_state_energy\n";

std::ofstream dr_file("convergence_dr.dat");

for(double test_dr : {0.5, 0.4, 0.3, 0.2, 0.15}) {
    matrix Htest = hydrogen_hamiltonian(rmax, test_dr);

    pp::Jacobi test_evd(Htest);

    vector ew = test_evd.w;
    matrix eV = test_evd.V;

    sort_eigenpairs(ew, eV);

    std::cout << test_dr << "      " << ew[0] << '\n';
    dr_file << test_dr << " " << ew[0] << '\n';
}

std::cout << "\nConvergence with fixed dr = " << dr << " and changing rmax:\n";
std::cout << "rmax      ground_state_energy\n";

std::ofstream rmax_file("convergence_rmax.dat");

for(double test_rmax : {4.0, 6.0, 8.0, 10.0, 12.0}) {
    matrix Htest = hydrogen_hamiltonian(test_rmax, dr);

    pp::Jacobi test_evd(Htest);

    vector ew = test_evd.w;
    matrix eV = test_evd.V;

    sort_eigenpairs(ew, eV);

    std::cout << test_rmax << "      " << ew[0] << '\n';
    rmax_file << test_rmax << " " << ew[0] << '\n';
}
    std::cout << "\nConvergence with fixed dr = " << dr << " and changing rmax:\n";
    std::cout << "rmax      ground_state_energy\n";

    for(double test_rmax : {4.0, 6.0, 8.0, 10.0, 12.0}) {
        matrix Htest = hydrogen_hamiltonian(test_rmax, dr);

        pp::Jacobi test_evd(Htest);

        vector ew = test_evd.w;
        matrix eV = test_evd.V;

        sort_eigenpairs(ew, eV);

        std::cout << test_rmax << "      " << ew[0] << '\n';
    }

    std::cout << '\n';
}

void part_C()
{
    std::cout << "C. Scaling check for Jacobi diagonalization\n";
    std::cout << "n      time_seconds\n";

    for(int n : {10, 15, 20, 25, 30}) {
        matrix A = random_symmetric_matrix(n);

        auto start = std::chrono::high_resolution_clock::now();

        pp::Jacobi evd(A);

        auto stop = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = stop - start;

        std::cout << n << "      " << elapsed.count() << '\n';
    }
}

int main(int argc, char** argv)
{
    std::cout << std::setprecision(12);

    Options opt = read_options(argc, argv);

    part_A();
    part_B(opt.rmax, opt.dr, opt.nstates);
    part_C();

    return 0;
}