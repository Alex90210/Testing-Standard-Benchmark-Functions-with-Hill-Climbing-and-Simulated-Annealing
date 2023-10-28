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

std::string next_neighbour(const std::string& binary_string, unsigned index) {

    if (index >= binary_string.length())
        index = get_random_unsigned(0, binary_string.length());
    std::string copy_string = binary_string;

    copy_string[index] = (copy_string[index] == '1') ? '0' : '1';
    ++index;

    return copy_string;
}

double hill_climbing(const double& interval_start, const double& interval_end, double epsilon,
                     unsigned number_of_dimensions, unsigned iterations, const std::string& mode,
                     double (*calculate_function)(const std::vector<double>& vec)) {

    double best_string_value_solution {10000};
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
        double solution = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, this_iteration_random_string));
        if (solution < best_string_value_solution)
            best_string_value_solution = solution;
    }
    return best_string_value_solution;
}

std::string generateNeighbor(const std::string& currentSolution, int numBitsToFlip) {
    std::string neighbor = currentSolution;
    std::random_device rd; // Initialize a random seed
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, currentSolution.size() - 1);

    for (int i = 0; i < numBitsToFlip; ++i) {
        int randomPosition = dis(gen);
        neighbor[randomPosition] = (neighbor[randomPosition] == '0') ? '1' : '0'; // Flip the bit
    }

    return neighbor;
}



double simulated_annealing(const double& interval_start, const double& interval_end, double epsilon,
                           unsigned number_of_dimensions, unsigned iterations, double temperature,
                           double (*calculate_function)(const std::vector<double>& vec)) {

    double temp_copy = temperature;
    double best {1000000};
    double cooling_rate {0.99};
    for(unsigned ii {0}; ii < iterations; ++ii) {
        temperature = temp_copy;
        std::string best_binary_string = generate_binary_string(interval_start, interval_end, epsilon,
                                                                number_of_dimensions);
        double best_value = calculate_function(
                decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, best_binary_string));

        while (temperature > 0.000001) {

            bool no_solution{false};
            unsigned iterator{0};
            while (!no_solution) {

                std::string random_neighbour1 = generateNeighbor(best_binary_string, 1);
                unsigned fate = get_random_unsigned(0, 1);
                /*if (fate == 1)
                    random_neighbour1 = random_neighbour(interval_start, interval_end, epsilon, best_binary_string);
                else
                    random_neighbour1 = first_improvement(interval_start, interval_end, epsilon, number_of_dimensions,
                                                          best_binary_string, best_value, calculate_function);
                ++iterator; // works better than nothing*/
                double random_neighbour_value = calculate_function(
                        decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions,
                                             random_neighbour1));

                if (random_neighbour_value < best_value) {
                    best_binary_string = random_neighbour1;
                    best_value = random_neighbour_value;
                } else {

                    double delta = random_neighbour_value - best_value;
                    /*double acceptance_probability = 1 / std::pow(10, delta / current_temperature);*/
                    double acceptance_probability = exp(-delta / temperature);
                    if (get_random_double(0, 1) < acceptance_probability) {
                        best_binary_string = random_neighbour1;
                        best_value = random_neighbour_value;
                    } else
                        no_solution = true;

                }
            }
            temperature *= cooling_rate;
        }
        if (best_value < best)
            best = best_value;
    }
    return best;
}

double simulated_annealing_adaptive(const double& interval_start, const double& interval_end, double epsilon,
                                    unsigned number_of_dimensions, unsigned iterations, double temperature,
                                    double (*calculate_function)(const std::vector<double>& vec)) {

    double temp_copy = temperature;
    double best {1000000};
    double cooling_rate {0.99};
    double initial_temperature = temperature; // Store the initial temperature

    for(unsigned ii {0}; ii < iterations; ++ii) {
        temperature = temp_copy;
        std::string best_binary_string = generate_binary_string(interval_start, interval_end, epsilon,
                                                                number_of_dimensions);
        double best_value = calculate_function(
                decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, best_binary_string));

        int accepted_worse_solutions = 0; // Counter for accepted worse solutions
        int total_worse_solutions = 0;    // Counter for total worse solutions
        int max_inner_iterations = 10000; // Add a maximum number of inner iterations

        while (temperature > 0.000001 && max_inner_iterations > 0) {
            bool no_solution{false};
            unsigned iterator{0};

            while (!no_solution) {
                std::string random_neighbour1 = generateNeighbor(best_binary_string, 1);
                double random_neighbour_value = calculate_function(
                        decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions,
                                             random_neighbour1));

                if (random_neighbour_value < best_value) {
                    best_binary_string = random_neighbour1;
                    best_value = random_neighbour_value;
                } else {

                    double delta = random_neighbour_value - best_value;
                    double acceptance_probability = exp(-delta / temperature);

                    if (get_random_double(0, 1) < acceptance_probability) {
                        best_binary_string = random_neighbour1;
                        best_value = random_neighbour_value;
                        accepted_worse_solutions++; // Increase counter for accepted worse solutions
                    } else
                        no_solution = true;
                }

                total_worse_solutions++; // Increase counter for total worse solutions
                max_inner_iterations--; // Decrease inner iteration counter
            }

            // Adapt the cooling rate based on the acceptance rate
            double acceptance_rate = (double)accepted_worse_solutions / total_worse_solutions;
            if (acceptance_rate < 0.2) {
                cooling_rate *= 0.9; // Reduce cooling rate
            } else if (acceptance_rate > 0.4) {
                cooling_rate *= 1.1; // Increase cooling rate
            }

            temperature *= cooling_rate;
        }

        if (best_value < best)
            best = best_value;
    }

    return best;
}

double simulated_annealing_with_decay(const double& interval_start, const double& interval_end, double epsilon,
                                      unsigned number_of_dimensions, unsigned iterations, double initial_temperature,
                                      double (*calculate_function)(const std::vector<double>& vec)) {

    double best {1000000};
    double cooling_rate {0.99};

    for(unsigned ii {0}; ii < iterations; ++ii) {
        double temperature = initial_temperature * exp(-0.01 * ii); // Exponential decay function

        std::string best_binary_string = generate_binary_string(interval_start, interval_end, epsilon,
                                                                number_of_dimensions);
        double best_value = calculate_function(
                decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, best_binary_string));

        int max_inner_iterations = 50;

        while (temperature > 0.000001 && max_inner_iterations > 0) {
            bool no_solution{false};
            unsigned iterator{0};

            while (!no_solution) {
                unsigned index{0};
                std::string random_neighbour1 = generateNeighbor(best_binary_string, 1);
                index++;
                double random_neighbour_value = calculate_function(
                        decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions,
                                             random_neighbour1));

                if (random_neighbour_value < best_value) {
                    best_binary_string = random_neighbour1;
                    best_value = random_neighbour_value;
                } else {
                    double delta = random_neighbour_value - best_value;
                    double acceptance_probability = 1 / std::pow(10, delta / temperature);
                    // double acceptance_probability = exp(-delta / temperature);

                    if (get_random_double(0, 1) < acceptance_probability) {
                        best_binary_string = random_neighbour1;
                        best_value = random_neighbour_value;
                    } else
                        no_solution = true;
                }

                max_inner_iterations--;
            }

            temperature *= cooling_rate;
        }

        if (best_value < best)
            best = best_value;
    }

    return best;
}

double simulated_annealing_with_linear_decay(const double& interval_start, const double& interval_end, double epsilon,
                                             unsigned number_of_dimensions, unsigned iterations, double initial_temperature,
                                             double (*calculate_function)(const std::vector<double>& vec)) {

    double best {1000000};
    double cooling_rate {0.99};
    double temperature = initial_temperature;

    for(unsigned ii {0}; ii < iterations; ++ii) {
        // Linear decay function
        temperature = initial_temperature - (ii * (initial_temperature / iterations));

        std::string best_binary_string = generate_binary_string(interval_start, interval_end, epsilon,
                                                                number_of_dimensions);
        double best_value = calculate_function(
                decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, best_binary_string));

        int max_inner_iterations = 10000;

        while (temperature > 0.000001 && max_inner_iterations > 0) {
            bool no_solution{false};
            unsigned iterator{0};

            while (!no_solution) {
                std::string random_neighbour1 = generateNeighbor(best_binary_string, 2);
                double random_neighbour_value = calculate_function(
                        decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions,
                                             random_neighbour1));

                if (random_neighbour_value < best_value) {
                    best_binary_string = random_neighbour1;
                    best_value = random_neighbour_value;
                } else {
                    double delta = random_neighbour_value - best_value;
                    double acceptance_probability = exp(-delta / temperature);

                    if (get_random_double(0, 1) < acceptance_probability) {
                        best_binary_string = random_neighbour1;
                        best_value = random_neighbour_value;
                    } else
                        no_solution = true;
                }

                max_inner_iterations--;
            }

            temperature *= cooling_rate;
        }

        if (best_value < best)
            best = best_value;
    }

    return best;
}
