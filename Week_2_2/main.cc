#include <iostream>
#include <vector>
#include "hello.h"
#include"vec.h"
int main() {
    hello();

    double a = 1.0;
    double b = a;

    if (a == b) {
        std::cout << "a == b\n";
    } else {
        std::cout << "a != b\n";
    }

    std::vector<double> v{1, 2, 3};

    for (std::size_t i = 0; i < v.size(); ++i) {
        std::cout << v[i] << ' ';
    }
    std::cout << '\n';

    for (auto vi : v) std::cout << vi << ' ';
    std::cout << '\n';

    for (double vi : v) std::cout << vi << ' ';
    std::cout << '\n';

    // Doesn't modify v (vi is a copy)
    for (auto vi : v) vi = 6;
    for (auto vi : v) std::cout << vi << ' ';
    std::cout << '\n';

    // Modifies v (vi is a reference)
    for (auto& vi : v) vi = 6;
    for (const auto& vi : v) std::cout << vi << ' ';
    std::cout << "\nNow comes the while loop\n";

    std::size_t i = 0;
    while (i < v.size()) {
        std::cout << "v[" << i << "]=" << v[i] << '\n';
        i+=1;
    }
std::cout << "\nNow comes the do loop\n";
    i = 0;
    do {
        std::cout << "v[" << i << "]=" << v[i] << '\n';
        i+=1;
    } while (i < v.size());

std::cout << "\nNow we use the class vec we made";
    std::cout << '\n';
    pp::vec alpha {1,2,3};
    alpha.x=6;
    std::cout <<alpha.x <<' '<< alpha.y <<' '<< alpha.z <<'\n';
    return 0;
}
