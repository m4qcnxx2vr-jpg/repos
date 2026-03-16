#include "lsfit.h"
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <fstream>

int main() {
    pp::vector t, y, dy;

    t.push_back(1);   y.push_back(117);  dy.push_back(6);
    t.push_back(2);   y.push_back(100);  dy.push_back(5);
    t.push_back(3);   y.push_back(88);   dy.push_back(4);
    t.push_back(4);   y.push_back(72);   dy.push_back(4);
    t.push_back(6);   y.push_back(53);   dy.push_back(4);
    t.push_back(9);   y.push_back(29.5); dy.push_back(3);
    t.push_back(10);  y.push_back(25.2); dy.push_back(3);
    t.push_back(13);  y.push_back(15.2); dy.push_back(2);
    t.push_back(15);  y.push_back(11.1); dy.push_back(2);

    std::cout << "Original data:\n";
    std::cout << "i   t   y   dy\n";
    for(int i = 0; i < t.size(); i++){
        std::cout << i << "   "
                  << t[i] << "   "
                  << y[i] << "   "
                  << dy[i] << "\n";
    }

    pp::vector Y, dY;
    for(int i = 0; i < y.size(); i++){
        Y.push_back(std::log(y[i]));
        dY.push_back(dy[i] / y[i]);
    }

    std::cout << "\nTransformed data:\n";
    std::cout << "i   t   ln(y)   dln(y)\n";
    for(int i = 0; i < t.size(); i++){
        std::cout << i << "   "
                  << t[i] << "   "
                  << Y[i] << "   "
                  << dY[i] << "\n";
    }

    std::vector<pp::Func> fs{
        [](double){ return 1.0; },
        [](double z){ return z; }
    };

    auto [c, cov] = pp::lsfit(fs, t, Y, dY);

    double ln_a = c[0];
    double lambda = -c[1];
    double a = std::exp(ln_a);

    double dln_a = std::sqrt(cov(0,0));
    double dlambda = std::sqrt(cov(1,1));

    double T_half = std::log(2.0)/lambda;
    double dT_half = std::log(2.0)/(lambda*lambda) * dlambda;

    std::cout << "\nFit coefficients for ln(y) = c0 + c1*t:\n";
    std::cout << "c0 = " << c[0] << " +/- " << dln_a << "\n";
    std::cout << "c1 = " << c[1] << " +/- " << std::sqrt(cov(1,1)) << "\n";

    std::cout << "\nCovariance matrix:\n";
    cov.print();

    std::cout << "\nPhysical parameters in y(t) = a*exp(-lambda*t):\n";
    std::cout << "a = " << a << "\n";
    std::cout << "lambda = " << lambda << " +/- " << dlambda << " day^-1\n";
    std::cout << "Half-life T1/2 = " << T_half << " +/- " << dT_half << " days\n";

    std::cout << "\nModern value for 224Ra half-life: about 3.63 days\n";
    if(std::fabs(T_half - 3.63) <= dT_half)
        std::cout << "The modern value agrees within the estimated uncertainty.\n";
    else
        std::cout << "The modern value does NOT agree within the estimated uncertainty.\n";

    std::ofstream data("decay_data.txt");
    for(int i = 0; i < t.size(); i++){
        data << t[i] << " "
             << y[i] << " "
             << dy[i] << "\n";
    }

    std::ofstream params("fit_params.gpi");
    params << "c0 = " << c[0] << "\n";
    params << "c1 = " << c[1] << "\n";
    params << "dc0 = " << dln_a << "\n";
    params << "dc1 = " << std::sqrt(cov(1,1)) << "\n";
    params << "a = " << a << "\n";
    params << "lambda = " << lambda << "\n";

    std::cout << "\nCoefficient uncertainties:\n";
    std::cout << "dc0 = " << dln_a << "\n";
    std::cout << "dc1 = " << std::sqrt(cov(1,1)) << "\n";
    return 0;
}