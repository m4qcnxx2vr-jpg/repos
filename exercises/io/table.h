#pragma once

#include <iosfwd>
#include <vector>

void print_table_header(std::ostream& output);

void print_number_row(std::ostream& output, double x);

void print_numbers(std::ostream& output, const std::vector<double>& numbers);

void read_numbers_from_stream(std::istream& input, std::ostream& output);