#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include "select_neighbour_strategy.hpp"
#include <string>
#include <vector>
#include <chrono>

double hill_climbing(const double& interval_start, const double& interval_end, double epsilon,
                     unsigned number_of_dimensions, unsigned iterations, const std::string& mode,
                     double (*calculate_function)(const std::vector<double>& vec));
unsigned get_random_unsigned(unsigned min, unsigned max);
double get_random_double(double min, double max);
std::string random_neighbour(const double& interval_start, const double& interval_end, double epsilon, const std::string& binary_string);
std::string random_neighbour_one_bit(const std::string& binary_string);
std::string next_neighbour(const std::string& binary_string, unsigned index);
std::string generateNeighbor(const std::string& currentSolution, int bits_to_flip);
double simulated_annealing(const double& interval_start, const double& interval_end, double epsilon,
                           unsigned number_of_dimensions, unsigned iterations, double temperature,
                           double (*calculate_function)(const std::vector<double>& vec));
double simulated_annealing_adaptive(const double& interval_start, const double& interval_end, double epsilon,
                                    unsigned number_of_dimensions, unsigned iterations, double temperature,
                                    double (*calculate_function)(const std::vector<double>& vec));
double simulated_annealing_with_decay(const double& interval_start, const double& interval_end, double epsilon,
                                      unsigned number_of_dimensions, unsigned iterations, double initial_temperature,
                                      double (*calculate_function)(const std::vector<double>& vec));
double simulated_annealing_with_linear_decay(const double& interval_start, const double& interval_end, double epsilon,
                                             unsigned number_of_dimensions, unsigned iterations, double initial_temperature,
                                             double (*calculate_function)(const std::vector<double>& vec));

#endif