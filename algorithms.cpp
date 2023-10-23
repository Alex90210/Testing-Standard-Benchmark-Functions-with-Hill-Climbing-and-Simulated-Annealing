#include "algorithms.hpp"

unsigned get_random_unsigned(unsigned min, unsigned max) {
    std::random_device rd;
    std::mt19937 eng(rd());

    std::uniform_int_distribution<> distribution(min, max);

    return distribution(eng);
}

double get_random_double(double min, double max) {
    std::random_device rd;
    std::mt19937 eng(rd());

    std::uniform_real_distribution<> distribution(min, max);

    return distribution(eng);
}

std::string random_neighbour(const double& interval_start, const double& interval_end, double epsilon, const std::string& binary_string) {

    // at max only one bit should be flipped for every dimension, I don't know how to enforce the interval when flipping multiple bits
    // what is writen above is not true, the decoding puts the value in the interval, but the solution implemented below might not be the worst idea
    // I should experiment with multiple flipped bits for each dimension though
    std::string copy_string = binary_string;
    unsigned dim_len = D_binary_length(interval_start, interval_end, epsilon);
    unsigned max_number_of_flipped_bits = ceil(std::log2(binary_string.length()));
    unsigned random_max = get_random_unsigned(1, max_number_of_flipped_bits);

    unsigned dimensional_counter {0};
    for (size_t i {0}; i < random_max; ++i) {
        unsigned random_index = get_random_unsigned(dimensional_counter, binary_string.length());
        copy_string[random_index] = (copy_string[random_index] == '1') ? '0' : '1';
        dimensional_counter += dim_len;
    }

    return copy_string;
}

std::string random_neighbour_one_bit(const std::string& binary_string) {

    std::string copy_string = binary_string;
    unsigned random_index = get_random_unsigned(0, binary_string.length());

    copy_string[random_index] = (copy_string[random_index] == '1') ? '0' : '1';

    return copy_string;
}

std::string next_neighbour(const std::string& binary_string) {

    static int index{0};
    if (index >= binary_string.length())
        index = 0;
    std::string copy_string = binary_string;

    copy_string[index] = (copy_string[index] == '1') ? '0' : '1';
    ++index;

    return copy_string;
}

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
        } else if (mode == "BI"){ // best improvement
            while(!is_local_minimum) {

                std::string new_string = best_improvement(interval_start, interval_end, epsilon, number_of_dimensions, this_iteration_random_string, string_value, calculate_function);
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

double simulated_annealing(const double& interval_start, const double& interval_end, double epsilon,
                           unsigned number_of_dimensions, unsigned iterations, double temperature,
                           double (*calculate_function)(const std::vector<double>& vec)) {

    // ????????????????????
    const double cooling_rate {0.001};
    double current_temperature {temperature};
    double best {1000000};
    unsigned ii {0};
    while (ii < iterations) {
        std::string best_binary_string = generate_binary_string(interval_start, interval_end, epsilon, number_of_dimensions);
        double best_value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, best_binary_string));
        while (current_temperature > 0.00001) {
            bool no_solution {false};
            unsigned index {0};
            while (!no_solution) {
                // std::string random_neighbour = first_improvement(interval_start, interval_end, epsilon, number_of_dimensions, best_binary_string, best_value, calculate_function);
                std::string random_neighbour = next_neighbour(best_binary_string);


                double random_neighbour_value = calculate_function(
                        decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions,
                                             random_neighbour));

                if (random_neighbour_value < best_value) {
                    best_binary_string = random_neighbour;
                    best_value = random_neighbour_value;
                } else {
                    double delta = random_neighbour_value - best_value;
                    // acc probability is faulty
                    /*double acceptance_probability = 1 / std::pow(10, delta / current_temperature);*/
                    double acceptance_probability = exp(-delta / temperature);
                    if (get_random_double(0, 1) < acceptance_probability) {
                        best_binary_string = random_neighbour;
                        best_value = random_neighbour_value;
                    } else
                        no_solution = true;
                }
            }
            current_temperature *= temperature / (1 + cooling_rate * ii);;
        }
        if (best_value < best)
            best = best_value;
        ++ii;
    }

    return best;
}
