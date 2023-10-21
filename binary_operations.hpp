#ifndef BINARY_OPERATIONS_HPP
#define BINARY_OPERATIONS_HPP

#include "math_functions.hpp"
#include <iostream>
#include <string>
#include <random>

unsigned D_binary_length(const double& interval_start, const double& interval_end, double epsilon);
unsigned binary_to_decimal(const std::string& binary_string, const size_t& string_start, const size_t& string_end);
std::string generate_binary_string(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions);
std::vector<double> decode_binary_string(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions, const std::string& binary_string);

#endif
