#include "algorithms.hpp"

double hill_climbing(const double& interval_start, const double& interval_end, double epsilon,
                     unsigned number_of_dimensions, unsigned iterations, const std::string& mode,
                     double (*calculate_function)(const std::vector<double>& vec)) {

    double best_string_value_solution {1000000000};
    for (size_t i {0}; i < iterations; ++i) {
        bool is_local_minimum {false};
        std::string this_iteration_random_string = generate_binary_string(interval_start, interval_end, epsilon, number_of_dimensions);
        double string_value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, this_iteration_random_string));

        if (mode == "WI") { // worst improvement
            while(!is_local_minimum) {

                std::string new_string = worst_improvement(interval_start, interval_end, epsilon, number_of_dimensions, this_iteration_random_string, string_value, calculate_function);
                double new_string_value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, new_string));

                if (new_string_value < string_value) {
                    this_iteration_random_string = new_string;
                    string_value = new_string_value;
                }
                else
                    is_local_minimum = true;
            }
        } else if (mode == "FI") { // first improvement
            while(!is_local_minimum) {

                std::string new_string = first_improvement(interval_start, interval_end, epsilon, number_of_dimensions, this_iteration_random_string, string_value, calculate_function);
                double new_string_value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, new_string));

                if (new_string_value < string_value) {
                    this_iteration_random_string = new_string;
                    string_value = new_string_value;
                }
                else
                    is_local_minimum = true;
            }
        } else { // best improvement
            while(!is_local_minimum) {

                std::string new_string = best_improvement(interval_start, interval_end, epsilon, number_of_dimensions, this_iteration_random_string, calculate_function);
                double new_string_value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, new_string));

                if (new_string_value < string_value) {
                    this_iteration_random_string = new_string;
                    string_value = new_string_value;
                }
                else
                    is_local_minimum = true;
            }
        }

        if (calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, this_iteration_random_string)) < best_string_value_solution)
            best_string_value_solution = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, this_iteration_random_string));
    }
    return best_string_value_solution;
}

unsigned get_random_number(unsigned min, unsigned max) {
    std::random_device rd;
    std::mt19937 eng(rd());

    std::uniform_int_distribution<> distribution(min, max);

    return distribution(eng);
}

std::string random_neighbour(const double& interval_start, const double& interval_end, double epsilon, const std::string& binary_string) {

    // at max only one bit should be flipped for every dimension, I don't know how to enforce the interval when flipping multiple bits
    // what is writen above is not true, the decoding puts the value in the interval, but the solution implemented below might not be the worst idea
    // I should experiment with multiple flipped bits for each dimension though
    std::string copy_string = binary_string;
    unsigned dim_len = D_binary_length(interval_start, interval_end, epsilon);
    unsigned max_number_of_flipped_bits = ceil(std::log2(binary_string.length()));
    unsigned random_max = get_random_number(1, max_number_of_flipped_bits);

    unsigned dimensional_counter {0};
    for (size_t i {0}; i < random_max; ++i) {
        unsigned random_index = get_random_number(dimensional_counter, binary_string.length());
        copy_string[random_index] = (copy_string[random_index] == '1') ? '0' : '1';
        dimensional_counter += dim_len;
    }

    return copy_string;
}

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