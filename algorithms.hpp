#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include "select_neighbour_strategy.hpp"

double hill_climbing(const double& interval_start, const double& interval_end, double epsilon,
                     unsigned number_of_dimensions, unsigned iterations, const std::string& mode,
                     double (*calculate_function)(const std::vector<double>& vec));
unsigned get_random_number(unsigned min, unsigned max);
std::string random_neighbour(const double& interval_start, const double& interval_end, double epsilon, const std::string& binary_string);

/*double simulated_annealing(const double& interval_start, const double& interval_end, double epsilon,
                           unsigned number_of_dimensions, unsigned iterations, const std::string& mode,
                           double (*calculate_function)(const std::vector<double>& vec)) {

    std::string random_string = generate_binary_string(interval_start, interval_end, epsilon, number_of_dimensions);
    double random_string_value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, random_string));

    do {
        do {

        } while ();
    } while();

}*/

#endif