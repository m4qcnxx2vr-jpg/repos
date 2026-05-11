#include "table.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
    std::vector<double> numbers;

    std::string input_filename;
    std::string output_filename;

    bool using_command_line_numbers = false;
    bool using_file_streams = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-n" && i + 1 < argc) {
            numbers.push_back(std::stod(argv[i + 1]));
            using_command_line_numbers = true;
            i++;
        } else if (arg == "--input" && i + 1 < argc) {
            input_filename = argv[i + 1];
            using_file_streams = true;
            i++;
        } else if (arg == "--output" && i + 1 < argc) {
            output_filename = argv[i + 1];
            using_file_streams = true;
            i++;
        }
    }

    if (using_file_streams) {
        if (input_filename.empty() || output_filename.empty()) {
            std::cerr << "Error: both --input and --output are required.\n";
            return EXIT_FAILURE;
        }

        std::ifstream input_file(input_filename);
        std::ofstream output_file(output_filename);

        if (!input_file.is_open()) {
            std::cerr << "Error: could not open input file: "
                      << input_filename << "\n";
            return EXIT_FAILURE;
        }

        if (!output_file.is_open()) {
            std::cerr << "Error: could not open output file: "
                      << output_filename << "\n";
            return EXIT_FAILURE;
        }

        read_numbers_from_stream(input_file, output_file);

        return EXIT_SUCCESS;
    }

    if (using_command_line_numbers) {
        print_numbers(std::cout, numbers);
        return EXIT_SUCCESS;
    }

    read_numbers_from_stream(std::cin, std::cout);

    return EXIT_SUCCESS;
}