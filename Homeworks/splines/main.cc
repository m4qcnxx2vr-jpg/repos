#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "linterp.h"
#include "qspline.h"
#include "fqinterp.h"

// ---------------- Debug function for Task B ----------------
void debug_qspline(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::string& name)
{
    qspline s(x,y);

    std::cout << "----- " << name << " -----\n";

    std::cout << "b: ";
    for(int i = 0; i < s.n-1; i++) std::cout << s.b[i] << " ";
    std::cout << "\n";

    std::cout << "c: ";
    for(int i = 0; i < s.n-1; i++) std::cout << s.c[i] << " ";
    std::cout << "\n\n";
}

// ---------------- Main ----------------
int main(){

    // ================= Task A =================
    {
        std::vector<double> x, y;

        for(double xi = 0; xi <= 9.0; xi += 0.5){
            x.push_back(xi);
            y.push_back(std::cos(xi));
        }

        std::ofstream data("linspline_data.txt");
        for(std::size_t i = 0; i < x.size(); i++){
            data << x[i] << " " << y[i] << "\n";
        }

        std::ofstream curve("linspline_curve.txt");
        for(double z = x.front(); z <= x.back(); z += 0.02){
            curve << z << " "
                  << linterp(x,y,z) << " "
                  << linterpInteg(x,y,z) << "\n";
        }

        std::cout << "===== Task A: Linear spline =====\n";

        for(double z : {0.25, 1.3, 4.7, 8.8}){
            std::cout << "z = " << z
                      << "  linterp = " << linterp(x,y,z)
                      << "  linterpInteg = " << linterpInteg(x,y,z)
                      << "\n";
        }

        std::cout << "\n";
    }

    // ================= Task B =================
    {
        std::cout << "===== Task B: Quadratic spline =====\n";

        // ---- Debug tests (VERY IMPORTANT for assignment) ----
        std::vector<double> xd{1,2,3,4,5};

        debug_qspline(xd, {1,1,1,1,1}, "y = 1");
        debug_qspline(xd, {1,2,3,4,5}, "y = x");
        debug_qspline(xd, {1,4,9,16,25}, "y = x^2");

        // ---- Actual spline example (cos) ----
        std::vector<double> x, y;

        for(double xi = 0; xi <= 9.0; xi += 0.5){
            x.push_back(xi);
            y.push_back(std::cos(xi));
        }

        qspline s(x,y);

        std::ofstream data("qspline_data.txt");
        for(std::size_t i = 0; i < x.size(); i++){
            data << x[i] << " " << y[i] << "\n";
        }

        std::ofstream curve("qspline_curve.txt");
        for(double z = x.front(); z <= x.back(); z += 0.02){
            curve << z << " "
                  << s.eval(z) << " "
                  << s.deriv(z) << " "
                  << s.integ(z) << "\n";
        }

        std::cout << "Sample values:\n";
        for(double z : {0.25, 1.3, 4.7, 8.8}){
            std::cout << "z = " << z
                      << "  eval = " << s.eval(z)
                      << "  deriv = " << s.deriv(z)
                      << "  integ = " << s.integ(z)
                      << "\n";
        }

        std::cout << "\n";
    }

    // ================= Task C =================
    {
        std::cout << "===== Task C: Functional quadratic spline =====\n";

        std::vector<double> x, y;
        for(double xi = 0; xi <= 9.0; xi += 0.5){
            x.push_back(xi);
            y.push_back(std::cos(xi));
        }

        auto f = make_qspline(x, y);

        std::ofstream curve("fqinterp_curve.txt");
        for(double z = x.front(); z <= x.back(); z += 0.02){
            curve << z << " " << f(z) << "\n";
        }

        for(double z : {0.25, 1.3, 4.7, 8.8}){
            std::cout << "z = " << z
                      << "  f(z) = " << f(z)
                      << "\n";
        }

        std::cout << "\n";
    }

    return 0;
}