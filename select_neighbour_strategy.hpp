#ifndef SELECT_NEIGHBOUR_STRATEGY_HPP
#define SELECT_NEIGHBOUR_STRATEGY_HPP

#include "binary_operations.hpp"
/*#include "math_functions.hpp"*/

std::string best_improvement(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions,
                             const std::string& binary_string, double (*calculate_function)(const std::vector<double>& vec));
std::string worst_improvement(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions,
                              const std::string& binary_string, double string_value,
                              double (*calculate_function)(const std::vector<double>& vec));
std::string first_improvement(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions,
                              const std::string& binary_string, double string_value,
                              double (*calculate_function)(const std::vector<double>& vec));

#endif