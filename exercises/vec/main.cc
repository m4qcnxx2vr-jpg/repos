#include "vec.h"

#include <iostream>
#include <cmath>

void test_vec(std::string name, const vec& calculated, const vec& expected)
{
    std::cout << name << '\n';
    std::cout << "  calculated = " << calculated << '\n';
    std::cout << "  expected   = " << expected << '\n';

    if (approx(calculated, expected)) {
        std::cout << "  result     = PASS\n\n";
    } else {
        std::cout << "  result     = FAIL\n\n";
    }
}

void test_double(std::string name, double calculated, double expected)
{
    std::cout << name << '\n';
    std::cout << "  calculated = " << calculated << '\n';
    std::cout << "  expected   = " << expected << '\n';

    if (approx(calculated, expected)) {
        std::cout << "  result     = PASS\n\n";
    } else {
        std::cout << "  result     = FAIL\n\n";
    }
}

int main()
{
    vec a(1.0, 2.0, 3.0);
    vec b(4.0, 5.0, 6.0);

    std::cout << "Testing vec implementation\n\n";

    std::cout << "a = " << a << '\n';
    std::cout << "b = " << b << "\n\n";

    std::cout << "Addition:\n";
    std::cout << "  a + b = { 1 + 4, 2 + 5, 3 + 6 }\n";
    test_vec("  test", a + b, vec(5.0, 7.0, 9.0));

    std::cout << "Subtraction:\n";
    std::cout << "  b - a = { 4 - 1, 5 - 2, 6 - 3 }\n";
    test_vec("  test", b - a, vec(3.0, 3.0, 3.0));

    std::cout << "Vector times scalar:\n";
    std::cout << "  a * 2 = { 1 * 2, 2 * 2, 3 * 2 }\n";
    test_vec("  test", a * 2.0, vec(2.0, 4.0, 6.0));

    std::cout << "Scalar times vector:\n";
    std::cout << "  2 * a = { 2 * 1, 2 * 2, 2 * 3 }\n";
    test_vec("  test", 2.0 * a, vec(2.0, 4.0, 6.0));

    std::cout << "Division by scalar:\n";
    std::cout << "  b / 2 = { 4 / 2, 5 / 2, 6 / 2 }\n";
    test_vec("  test", b / 2.0, vec(2.0, 2.5, 3.0));

    std::cout << "Dot product:\n";
    std::cout << "  a dot b = 1 * 4 + 2 * 5 + 3 * 6\n";
    std::cout << "          = 4 + 10 + 18\n";
    test_double("  test", a.dot(b), 32.0);

    std::cout << "Cross product:\n";
    std::cout << "  a cross b = {\n";
    std::cout << "      2 * 6 - 3 * 5,\n";
    std::cout << "      3 * 4 - 1 * 6,\n";
    std::cout << "      1 * 5 - 2 * 4\n";
    std::cout << "  }\n";
    std::cout << "  a cross b = { 12 - 15, 12 - 6, 5 - 8 }\n";
    test_vec("  test", a.cross(b), vec(-3.0, 6.0, -3.0));

    std::cout << "Norm:\n";
    std::cout << "  norm(a) = sqrt(a dot a)\n";
    std::cout << "          = sqrt(1 * 1 + 2 * 2 + 3 * 3)\n";
    std::cout << "          = sqrt(1 + 4 + 9)\n";
    std::cout << "          = sqrt(14)\n";
    test_double("  test", a.norm(), std::sqrt(14.0));

    std::cout << "Default constructor:\n";
    std::cout << "  vec zero;\n";
    std::cout << "  expected zero vector = { 0, 0, 0 }\n";
    vec zero;
    test_vec("  test", zero, vec(0.0, 0.0, 0.0));

    std::cout << "Copy constructor:\n";
    std::cout << "  vec copy = a;\n";
    std::cout << "  expected copy = a = { 1, 2, 3 }\n";
    vec copy = a;
    test_vec("  test", copy, a);

    std::cout << "Copy assignment:\n";
    std::cout << "  assigned = b;\n";
    std::cout << "  expected assigned = b = { 4, 5, 6 }\n";
    vec assigned;
    assigned = b;
    test_vec("  test", assigned, b);

    std::cout << "Print method:\n";
    a.print("  a printed with print method: ");

    return 0;
}