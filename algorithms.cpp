#include "algorithms.hpp"

unsigned get_random_unsigned(unsigned min, unsigned max) {
    static std::random_device rd;
    static std::mt19937_64 eng(rd());

    std::uniform_int_distribution<unsigned> distribution(min, max);

    return distribution(eng);
}

double get_random_double(double min, double max) {
    static std::random_device rd;
    static std::mt19937_64 eng(rd());

    std::uniform_real_distribution<double> distribution(min, max);

    return distribution(eng);
}

std::string random_neighbour(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions, std::string& binary_string) {

    std::string copy_string = binary_string;
    unsigned dim_len = D_binary_length(interval_start, interval_end, epsilon);

    for (size_t i {0}; i < number_of_dimensions; ++i) {
        for (size_t j {i * dim_len}; j < (i + 1) * dim_len; j++) {
            unsigned index = get_random_unsigned(i * dim_len, (i + 1) * dim_len);
            copy_string[index] = (copy_string[index] == '1') ? '0' : '1';
        }
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

        if (mode == "UP") { // up best improvement
            std::string this_iteration_random_string_up = generate_binary_string_h1p(interval_start, interval_end, epsilon, number_of_dimensions);
            while(!is_local_minimum) {

                std::string new_string = best_improvement(interval_start, interval_end, epsilon, number_of_dimensions, this_iteration_random_string, string_value, calculate_function);
                double new_string_value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, new_string));

                if (new_string_value > string_value) {
                    this_iteration_random_string_up = new_string;
                    string_value = new_string_value;
                }
                else
                    is_local_minimum = true;
            }
        }
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

std::string generate_neighbor_n_flipped_bits(const std::string& currentSolution, int numBitsToFlip) {
    std::string neighbor = currentSolution;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, currentSolution.size() - 1);

    for (int i = 0; i < numBitsToFlip; ++i) {
        int randomPosition = dis(gen);
        neighbor[randomPosition] = (neighbor[randomPosition] == '0') ? '1' : '0';
    }

    return neighbor;
}

double simulated_annealing(const double& interval_start, const double& interval_end, double epsilon,
                           unsigned number_of_dimensions, unsigned iterations, double temperature,
                           double (*calculate_function)(const std::vector<double>& vec)) {

    double best {1000000};
    double cooling_rate {0.999};

    for (unsigned ii {0}; ii < iterations; ++ii) {
        //double c_temperature = temperature;
        double c_temperature = temperature - (ii * (temperature / iterations));
        //double temperature = initial_temperature * exp(-0.01 * ii);
        std::string best_binary_string = generate_binary_string(interval_start, interval_end, epsilon,
                                                                number_of_dimensions);
        double best_value = calculate_function(
                decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, best_binary_string));

        int max_inner_iterations {200000};
        while (c_temperature >= 0.000000001 && max_inner_iterations > 0) {

            int no_solution{0};
            while (no_solution < 5) {
                std::string random_neighbour1 = generate_neighbor_n_flipped_bits(best_binary_string, 5);

                //std::string random_neighbour1 = random_neighbour(interval_start, interval_end, epsilon, number_of_dimensions, best_binary_string);
                double random_neighbour_value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions,random_neighbour1));

                if (random_neighbour_value < best_value) {
                    best_binary_string = random_neighbour1;
                    best_value = random_neighbour_value;
                } else {

                    double delta = random_neighbour_value - best_value;
                    double acceptance_probability = exp(-delta / c_temperature);
                    //double acceptance_probability = 1 / std::pow(10, delta / temperature);

                    if (get_random_double(0, 1) < acceptance_probability) {
                        best_binary_string = random_neighbour1;
                        best_value = random_neighbour_value;
                    } else
                        ++no_solution;
                }
                max_inner_iterations--;
            }
            c_temperature *= cooling_rate;
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

    for(unsigned ii {0}; ii < iterations; ++ii) {
        temperature = temp_copy;
        std::string best_binary_string = generate_binary_string(interval_start, interval_end, epsilon,
                                                                number_of_dimensions);
        double best_value = calculate_function(
                decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, best_binary_string));

        int accepted_worse_solutions = 0;
        int total_worse_solutions = 0;
        int max_inner_iterations = 10000;

        while (temperature > 0.000001 && max_inner_iterations > 0) {
            bool no_solution{false};

            while (!no_solution) {
                std::string random_neighbour1 = generate_neighbor_n_flipped_bits(best_binary_string, 1);
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

            double acceptance_rate = (double)accepted_worse_solutions / total_worse_solutions;
            if (acceptance_rate < 0.2) {
                cooling_rate *= 0.9;
            } else if (acceptance_rate > 0.4) {
                cooling_rate *= 1.1;
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
    double cooling_rate {0.9955};

    for(unsigned ii {0}; ii < iterations; ++ii) {
        double temperature = initial_temperature * exp(-0.01 * ii);

        std::string best_binary_string = generate_binary_string(interval_start, interval_end, epsilon,
                                                                number_of_dimensions);
        double best_value = calculate_function(
                decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, best_binary_string));

        int max_inner_iterations = 10000;

        while (temperature > 0.00001 && max_inner_iterations > 0) {
            bool no_solution{false};

            while (!no_solution) {
                unsigned index{0};
                std::string random_neighbour1 = generate_neighbor_n_flipped_bits(best_binary_string, 2);
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
    double cooling_rate {0.999};

    for (unsigned ii {0}; ii < iterations; ++ii) {
        double temperature = initial_temperature;
        //double temperature = initial_temperature - (ii * (initial_temperature / iterations));
        //double temperature = initial_temperature * exp(-0.01 * ii);
        std::string best_binary_string = generate_binary_string(interval_start, interval_end, epsilon,
                                                                number_of_dimensions);
        double best_value = calculate_function(
                decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, best_binary_string));
        int max_inner_iterations {200000};

        while (temperature > 0.00000001 && max_inner_iterations > 0) {

            int no_solution{0};
            while (no_solution < 5) {
                std::string random_neighbour1 = generate_neighbor_n_flipped_bits(best_binary_string, 5);

                //std::string random_neighbour1 = random_neighbour(interval_start, interval_end, epsilon, number_of_dimensions, best_binary_string);
                double random_neighbour_value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions,random_neighbour1));

                if (random_neighbour_value < best_value) {
                    best_binary_string = random_neighbour1;
                    best_value = random_neighbour_value;
                } else {

                    double delta = random_neighbour_value - best_value;
                    double acceptance_probability = exp(-delta / temperature);
                    //double acceptance_probability = 1 / std::pow(10, delta / temperature);

                    if (get_random_double(0, 1) < acceptance_probability) {
                        best_binary_string = random_neighbour1;
                        best_value = random_neighbour_value;
                    } else
                        ++no_solution;
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
