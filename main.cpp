#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <chrono>

// hill climbing algorithm steps:
// 1. generate random bitstring
// 2. decode bitstring and evaluate it
// 3. do-while improve strategy (the best improvement, the least improvement, first improvement) on generated neighbourhood

unsigned D_binary_length(const double& interval_start, const double& interval_end, double epsilon) {
    unsigned pow_epsilon = 1 / epsilon;
    unsigned dim_number_of_bits = std::ceil(std::log2((interval_end - interval_start) * pow_epsilon));
    return dim_number_of_bits;
}

std::string generate_binary_string(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);

    unsigned pow_epsilon = 1 / epsilon;
    unsigned dim_number_of_bits = std::ceil(std::log2((interval_end - interval_start) * pow_epsilon));
    unsigned number_of_bits = number_of_dimensions * dim_number_of_bits;

    std::string generated_string;
    for (size_t i = 0; i < number_of_bits; ++i) {
        generated_string += (dis(gen) == 0) ? '1' : '0';
    }

    return generated_string;
}

unsigned binary_to_decimal(const std::string& binary_string, const size_t& string_start, const size_t& string_end) {

    unsigned decimal_value {0};
    for (size_t i {string_start}; i < string_end; ++i) {
        decimal_value *= 2;
        decimal_value += binary_string[i] - '0';
    }

    return decimal_value;
}

std::vector<double> decode_binary_string(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions, const std::string& binary_string) {

    // x will be between 0 and 2^n - 1, n is the length of the binary string
    std::vector <double> dimensional_values;
    unsigned dim_length = D_binary_length(interval_start, interval_end, epsilon);

    for (size_t i {0}; i < number_of_dimensions; ++i) {

        unsigned xb_value = binary_to_decimal(binary_string, dim_length * i, dim_length * (i + 1));
        double x_value = xb_value / (pow(2, dim_length) - 1);

        x_value *= (interval_end - interval_start);
        x_value += interval_start; // x is now between interval_start and interval_end
        dimensional_values.push_back(x_value);
    }

    return dimensional_values;
}

double dejong1_function(const std::vector<double>& vec) {
    double sum {0.0};
    for (auto x : vec)
        sum += x * x;
    return sum;
}

double schwefels_function(const std::vector<double>& vec) {
    double sum = 0.0;
    for (auto x : vec) {
        sum += -x * sin(sqrt(fabs(x)));
    }
    return sum;
}

double rastrigins_function(const std::vector<double>& vec) {
    double sum = 0.0;
    for (auto x : vec) {
        sum += x * x - 10 * cos(2 * M_PI * x);
    }
    return 10 * vec.size() + sum;
}

double michalewiczs_function(const std::vector<double>& vec) {
    double sum = 0.0;
    for (unsigned i = 0; i < vec.size(); ++i) {
        sum += sin(vec[i]) * pow(sin((i + 1) * vec[i] * vec[i] / M_PI), 20.0);
    }
    return -sum;
}

std::string best_improvement(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions,
                             const std::string& binary_string, double (*calculate_function)(const std::vector<double>& vec)) {

    int index {-1};
    double best_value {std::numeric_limits<double>::max()};
    std::string copy_string = binary_string;
    for (int i {0}; i < copy_string.length(); ++i) {
        copy_string[i] = (copy_string[i] == '1') ? '0' : '1';
        double value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, copy_string));

        if (value < best_value) {
            best_value = value;
            index = i;
        }
        copy_string[i] = (copy_string[i] == '1') ? '0' : '1';
    }

    if (best_value < calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, copy_string)))
        copy_string[index] = (copy_string[index] == '1') ? '0' : '1';

    return copy_string;
}

std::string worst_improvement(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions,
                             const std::string& binary_string, double string_value,
                             double (*calculate_function)(const std::vector<double>& vec)) {

    // something needs to be fixed here, the function executes way to fast
    int index {-1};
    double current_worst_but_better {-1000000000};
    int dimensional_length = D_binary_length(interval_start, interval_end, epsilon);
    std::string copy_string = binary_string;
    for (int i {0}; i < number_of_dimensions; ++i) { // this has to skip the least significant bits
        for (int j {dimensional_length * i}; j < (dimensional_length * (i + 1)) / 2; ++j) {
            copy_string[i] = (copy_string[i] == '1') ? '0' : '1';
            double value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, copy_string));

            if (value < string_value && value > current_worst_but_better) {
                current_worst_but_better = value;
                index = i;
            }
            copy_string[i] = (copy_string[i] == '1') ? '0' : '1';
        }
    }

    if (current_worst_but_better < calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, copy_string)))
        copy_string[index] = (copy_string[index] == '1') ? '0' : '1';

    return copy_string;
}

std::string first_improvement(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions,
                             const std::string& binary_string, double string_value,
                             double (*calculate_function)(const std::vector<double>& vec)) {

    double best_value {string_value};
    std::string copy_string = binary_string;
    for (int i {0}; i < copy_string.length(); ++i) {
        copy_string[i] = (copy_string[i] == '1') ? '0' : '1';
        double value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, copy_string));

        if (value < best_value) {
            best_value = value;
            break;
        }
        copy_string[i] = (copy_string[i] == '1') ? '0' : '1';
    }

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

int main () {

    auto start = std::chrono::high_resolution_clock::now();

    // there is something wrong in the worst improvement cased

    std::string mode{"WI"};
    unsigned number_of_dimensions {2};
    double epsilon {0.01};
    unsigned iterations {1000};
    double temperature {1000};

    double interval_start {-5.12};
    double interval_end {5.12};
    std::string test = generate_binary_string(interval_start, interval_end, epsilon, number_of_dimensions);
    std::string neighbour = random_neighbour(interval_start, interval_end,epsilon, test);
    std::cout << test << std::endl;
    std::cout << neighbour << std::endl;

    /*// De Jong 1
    // must be 0 (for 30 dimensions)

    double interval_start {-5.12};
    double interval_end {5.12};
    double best_d = hill_climbing(interval_start, interval_end, epsilon, number_of_dimensions, iterations, mode, dejong1_function);
    std::cout << std::fixed << std::setprecision(5) << "De Jung 1: " << best_d << std::endl;

    // Schwefel
    // must be under -10000 (for 30 dimensions)

    interval_start = -500;
    interval_end = 500;
    double best_s = hill_climbing(interval_start, interval_end, epsilon, number_of_dimensions, iterations, mode, schwefels_function);
    std::cout << std::fixed << std::setprecision(5) << "Schwefel: " << best_s << std::endl;

    // Rastrigin
    // must be under -25 (for 30 dimensions)

    interval_start = -5.12;
    interval_end = 5.12;
    double best_r = hill_climbing(interval_start, interval_end, epsilon, number_of_dimensions, iterations, mode, rastrigins_function);
    std::cout << std::fixed << std::setprecision(5) << "Rastrigin: " << best_r << std::endl;

    // Michaleiwcz
    // must be under -25 (for 30 dimensions)

    interval_start = 0;
    interval_end = M_PI;
    double best_m = hill_climbing(interval_start, interval_end, epsilon, number_of_dimensions, iterations, mode, michalewiczs_function);
    std::cout << std::fixed << std::setprecision(5) << "Michalewicz: " << best_m << std::endl;*/

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Program executed in: " << duration.count() << " seconds." << std::endl;

    return 0;
}

/*double hill_climbing_variant1() {
    // I could find the optimum for every dimension, this could be better if there is no relationship between these dimensional values
    for (size_t i {0}; i < number_of_dimensions; ++i) {
        for (size_t j = size_of_dim * i; j < size_of_dim * (i + 1); ++j) {
            unsigned string_to_decimal = binary_to_decimal(generated_string, size_of_dim * i, size_of_dim * (i + 1));
            // double current_value = decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, )
            double dimensional_value = de_jong_1()
            if (generated_string[j] == '0')
                generated_string[j] = '1';
            else generated_string[j] = '0';
            // .............
        }
    }
}*/

/*std::vector<double> generate_neighbourhood_initial_attempt(double interval_start, double interval_end, double epsilon, unsigned number_of_dimensions, const std::string& binary_string) {
    std::string flipped_binary_string = flip_bits(binary_string);
    std::vector<double> dim_values = decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, binary_string);
    // the de jung 1 value is needed for calculating the neighbours
    // double function_value = de_jong_1(dim_values); // affects the modularity of the function, this needs to be changed
    // de jung 1 needs to be calculated for every dimension
    unsigned interval_partitioner {2};
    double interval_len = interval_end - interval_start;
    std::vector<std::vector<double>> neighbours_in_all_dimensions;

    // 10111-00111-10000-10101

    unsigned size_of_dim = D_binary_length(interval_start, interval_end, epsilon);
    // std::cout << std::endl << flipped_binary_string.length() << std::endl;
    std::vector<double> this_dim_vector;
    for (size_t i {0}; i < number_of_dimensions; ++i) {
        for (size_t j = size_of_dim * i; j < size_of_dim * (i + 1); ++j) {
            if(flipped_binary_string[j] == '1') {
                this_dim_vector.push_back(dim_values[i] + interval_len / interval_partitioner);
                interval_partitioner *= 2;
            } else {
                this_dim_vector.push_back(dim_values[i] - interval_len / interval_partitioner);
                interval_partitioner *= 2;
            }
        }
        neighbours_in_all_dimensions.push_back(this_dim_vector);
        this_dim_vector.clear();
        interval_partitioner = 2;
    }

    std::vector<double> temporary_vec;
    std::vector<double> neighbourhood_functional_values;
    for (size_t j {0}; j < size_of_dim; ++j) {
        for (size_t i {0}; i < number_of_dimensions; ++i) {
            temporary_vec.push_back(neighbours_in_all_dimensions[i][j]);
        }
        neighbourhood_functional_values.push_back(de_jong_1(temporary_vec));
        temporary_vec.clear();
    }
    return neighbourhood_functional_values;
}*/

/*std::string flip_bits(const std::string& binary_string) { // this goes very well with the useless first generate_neighbourhood
    // function I wrote
    std::string flipped_string;
    for(auto i: binary_string) {
        if (i == '1')
            flipped_string += '0';
        else
            flipped_string += '1';
    }
    return flipped_string;
}*/

/*double best_improvement(const std::vector<double>& vec, const double& best_solution) {
    double solution {best_solution};
    for (auto i : vec)
        if(i < solution)
            solution = i;

    return solution;
}*/