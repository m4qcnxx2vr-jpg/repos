#include "solver.h"
#include <iostream>
#include <fstream>
#include <numbers>
#include <cmath>
#include <vector>

using namespace pp;

double epsilon = 0.0;

// ==========================================
// Planetary ODE
// ==========================================
vector orbit(double phi, const vector& y) {
    (void)phi; 
    vector dydphi(2);
    dydphi[0] = y[1];
    dydphi[1] = 1.0 - y[0] + epsilon * y[0] * y[0];
    return dydphi;
}

// ==========================================
// 3-Body Newtonian ODE
// ==========================================
vector three_body(double t, const vector& z) {
    (void)t;
    vector dzdt(12);

    for(int i = 0; i < 6; ++i) dzdt[6 + i] = z[i];

    for(int i = 0; i < 3; ++i) {
        double ax = 0, ay = 0;
        double xi = z[6 + 2*i], yi = z[7 + 2*i];

        for(int j = 0; j < 3; ++j) {
            if(i == j) continue;
            double dx = z[6 + 2*j] - xi;
            double dy = z[7 + 2*j] - yi;
            double distSq = dx*dx + dy*dy;
            double distInv3 = 1.0 / (distSq * std::sqrt(distSq));
            
            ax += dx * distInv3;
            ay += dy * distInv3;
        }
        dzdt[2*i]     = ax;
        dzdt[2*i + 1] = ay;
    }
    return dzdt;
}

void run_orbit(const std::string& file, double eps_val, double u0, double u0p, double tol) {
    std::ofstream out(file);
    epsilon = eps_val;
    
    // Fixed initialization
    vector y0(2);
    y0[0] = u0;
    y0[1] = u0p;

    auto [phi, y] = driver(orbit, 0.0, 20 * std::numbers::pi, y0, 0.01, tol, tol);
    for (size_t i = 0; i < phi.size(); i++) {
        if (y[i][0] > 1e-9) {
            double r = 1.0 / y[i][0];
            out << phi[i] << " " << r * std::cos(phi[i]) << " " << r * std::sin(phi[i]) << "\n";
        }
    }
    std::cout << "Created " << file << " with " << phi.size() << " points.\n";
}

void run_figure8(const std::string& file) {
    std::ofstream out(file);
    
    double x1 = 0.97000436, y1 = -0.24308753;
    double v1x = 0.466203685, v1y = 0.43236573;

    // Fixed initialization
    vector z0(12);
    z0[0] = v1x;    z0[1] = v1y;    // v1
    z0[2] = v1x;    z0[3] = v1y;    // v2
    z0[4] = -2*v1x; z0[5] = -2*v1y; // v3
    z0[6] = x1;     z0[7] = y1;     // r1
    z0[8] = -x1;    z0[9] = -y1;    // r2
    z0[10] = 0;     z0[11] = 0;     // r3

    auto [t, z] = driver(three_body, 0.0, 6.3259, z0, 0.1, 1e-5, 1e-5);

    for (size_t i = 0; i < t.size(); i++) {
        out << t[i] << " " << z[i][6] << " " << z[i][7] << " " 
            << z[i][8] << " " << z[i][9] << " " 
            << z[i][10] << " " << z[i][11] << "\n";
    }
    std::cout << "Created " << file << " (Fig-8) with " << t.size() << " points.\n";
}

int main() {
    run_orbit("1.txt", 0.0, 1.000001, 0.0, 1e-9);
    run_orbit("2.txt", 0.0, 1.0, -0.5, 1e-5);
    run_orbit("3.txt", 0.01, 1.0, -0.5, 1e-5);
    run_figure8("8.txt");
    return 0;
}