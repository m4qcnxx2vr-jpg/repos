#include "matrix.h"
#include "jacobi.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <cmath>

// -------- helpers --------

static void die_usage(){
    std::cerr << "Usage: ./main -rmax <number> -dr <number> [-nstates <int>]\n";
    std::cerr << "Example: ./main -rmax 10 -dr 0.3\n";
    std::exit(1);
}

static bool read_arg(int argc, char** argv, const std::string& key, double& val){
    for(int i=1;i<argc-1;i++){
        if(std::string(argv[i])==key){
            val = std::stod(argv[i+1]);
            return true;
        }
    }
    return false;
}

static bool read_arg_int(int argc, char** argv, const std::string& key, int& val){
    for(int i=1;i<argc-1;i++){
        if(std::string(argv[i])==key){
            val = std::stoi(argv[i+1]);
            return true;
        }
    }
    return false;
}

static double exact_energy(int n){
    return -0.5 / (double)(n*n);
}

static void sort_eigenpairs(pp::vector& w, pp::matrix& V){
    int n = w.size();
    std::vector<int> idx(n);
    for(int i=0;i<n;i++) idx[i]=i;

    std::sort(idx.begin(), idx.end(), [&](int a, int b){ return w[a] < w[b]; });

    pp::vector ws(n);
    pp::matrix Vs(n,n,0.0);
    for(int k=0;k<n;k++){
        ws[k] = w[idx[k]];
        for(int i=0;i<n;i++) Vs(i,k) = V(i, idx[k]);
    }
    w = ws;
    V = Vs;
}

static void normalize_on_grid(pp::vector& f, double dr){
    // normalize so sum_i |f_i|^2 dr = 1
    double s = 0.0;
    for(int i=0;i<f.size();i++) s += f[i]*f[i];
    s = std::sqrt(s * dr);
    if(s==0.0) return;
    for(int i=0;i<f.size();i++) f[i] /= s;
}

// -------- main --------

int main(int argc, char** argv){
    double rmax = 0.0;
    double dr   = 0.0;
    int nstates = 4; // how many lowest states to print/save

    if(!read_arg(argc, argv, "-rmax", rmax)) die_usage();
    if(!read_arg(argc, argv, "-dr", dr))     die_usage();
    read_arg_int(argc, argv, "-nstates", nstates);

    if(rmax <= 0 || dr <= 0) die_usage();

    int npoints = (int)(rmax/dr) - 1;
    if(npoints < 2){
        std::cerr << "Error: npoints too small. Increase rmax or decrease dr.\n";
        return 1;
    }

    std::cout << "rmax=" << rmax << "  dr=" << dr << "  npoints=" << npoints << "\n";

    // grid
    pp::vector r(npoints);
    for(int i=0;i<npoints;i++) r[i] = dr*(i+1);

    // Hamiltonian
    pp::matrix H(npoints, npoints, 0.0);

    // your snippet uses: (-0.5/dr/dr) times the tridiagonal with (-2 on diag, 1 offdiag)
    double t = -0.5/(dr*dr);

    for(int i=0;i<npoints-1;i++){
        H(i,i)     = -2*t;
        H(i,i+1)   =  1*t;
        H(i+1,i)   =  1*t;
    }
    H(npoints-1, npoints-1) = -2*t;

    for(int i=0;i<npoints;i++){
        H(i,i) += -1.0/r[i];
    }

    // Diagonalize
    pp::Jacobi evd(H, 1e-12);

    // Sort eigenpairs
    sort_eigenpairs(evd.w, evd.V);

    // Print lowest bound energies (negative ones)
    std::cout << "\nLowest eigenvalues (Hartree) and exact En=-1/(2n^2):\n";
    std::cout << std::setw(6) << "n"
              << std::setw(18) << "numeric"
              << std::setw(18) << "exact"
              << std::setw(18) << "abs error"
              << "\n";

    int printed = 0;
    for(int k=0;k<npoints && printed<nstates;k++){
        double ek = evd.w[k];
        if(ek >= 0) break; // bound states are negative
        int n = printed + 1;
        double eex = exact_energy(n);
        std::cout << std::setw(6) << n
                  << std::setw(18) << std::setprecision(10) << ek
                  << std::setw(18) << std::setprecision(10) << eex
                  << std::setw(18) << std::setprecision(10) << std::fabs(ek - eex)
                  << "\n";
        printed++;
    }

    if(printed==0){
        std::cout << "\nNo negative eigenvalues found. Try larger rmax and/or smaller dr.\n";
    }

    // Save eigenfunctions (reduced radial f(r)) for first few bound states
    // Each file: state_n.dat with columns: r  f(r)
    int saved = 0;
    for(int k=0;k<npoints && saved<nstates;k++){
        if(evd.w[k] >= 0) break;

        pp::vector f(npoints);
        for(int i=0;i<npoints;i++) f[i] = evd.V(i,k); // kth eigenvector column
        normalize_on_grid(f, dr);

        int n = saved + 1;
        std::string fname = "state_" + std::to_string(n) + ".dat";
        std::ofstream out(fname);
        for(int i=0;i<npoints;i++){
            out << r[i] << " " << f[i] << "\n";
        }
        out.close();

        saved++;
    }

    std::cout << "\nWrote eigenfunctions to files: state_1.dat, state_2.dat, ... (up to nstates)\n";
    std::cout << "Each file has columns: r  f(r)\n";

    return 0;
}