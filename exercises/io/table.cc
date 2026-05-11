#include "table.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

void print_table_header(std::ostream& output)
{
    output << std::setw(15) << "x"
           << std::setw(15) << "sin(x)"
           << std::setw(15) << "cos(x)"
           << "\n";

    output << std::setw(15) << "-------------"
           << std::setw(15) << "-------------"
           << std::setw(15) << "-------------"
           << "\n";
}

void print_number_row(std::ostream& output, double x)
{
    output << std::setw(15) << std::scientific << std::setprecision(6) << x
           << std::setw(15) << std::scientific << std::setprecision(6) << std::sin(x)
           << std::setw(15) << std::scientific << std::setprecision(6) << std::cos(x)
           << "\n";
}

void print_numbers(std::ostream& output, const std::vector<double>& numbers)
{
    print_table_header(output);

    for (double x : numbers) {
        print_number_row(output, x);
    }
}

void read_numbers_from_stream(std::istream& input, std::ostream& output)
{
    print_table_header(output);

    double x;

    while (input >> x) {
        print_number_row(output, x);
    }
}