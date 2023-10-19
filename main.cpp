#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

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

double de_jong_1(const std::vector<double>& vec) {

    double sum {0.0};
    for (auto x : vec)
        sum += x * x;
    return sum;
}

double best_improvement(const std::vector<double>& vec, const double& best_solution) {

    double solution {best_solution};
    for (auto i : vec)
        if(i < solution)
            solution = i;

    return solution;
}

std::string generate_neighbourhood_and_select(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions, const std::string& binary_string) {

    int index {-1};
    double best_value {1000000000};
    std::string copy_string = binary_string;
    for (int i {0}; i < copy_string.length(); ++i) {
        copy_string[i] = (copy_string[i] == '1') ? '0' : '1';
        double value = de_jong_1(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, copy_string));

        if (value < best_value) {
            best_value = value;
            index = i;
        }
        copy_string[i] = (copy_string[i] == '1') ? '0' : '1';
    }

    if (best_value < de_jong_1(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, copy_string)))
        copy_string[index] = (copy_string[index] == '1') ? '0' : '1';

    return copy_string;
}

double hill_climbing(const double& interval_start, const double& interval_end, double epsilon,
                     unsigned number_of_dimensions, unsigned iterations) {

    double best_string_value_solution {1000000000};
    for (size_t i {0}; i < iterations; ++i) {
        bool is_local_minimum {false};
        std::string this_iteration_random_string = generate_binary_string(interval_start, interval_end, epsilon, number_of_dimensions);
        double string_value = de_jong_1(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, this_iteration_random_string));

        while(!is_local_minimum) {
            std::string new_string = generate_neighbourhood_and_select(interval_start, interval_end, epsilon, number_of_dimensions, this_iteration_random_string);
            double new_string_value = de_jong_1(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, new_string));

            if (new_string_value < string_value) {
                this_iteration_random_string = new_string;
                string_value = new_string_value;
            }
            else
                is_local_minimum = true;

        }

        if (de_jong_1(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, this_iteration_random_string)) < best_string_value_solution)
            best_string_value_solution = de_jong_1(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, this_iteration_random_string));
    }
    return best_string_value_solution;
}

int main () {

    double interval_start = -5.12;
    double interval_end = 5.12;
    double epsilon = 0.001;
    unsigned number_of_dimensions = 100;
    unsigned iterations {1};

    double best = hill_climbing(interval_start, interval_end, epsilon, number_of_dimensions, iterations);
    std::cout << std::fixed << std::setprecision(5) << best;

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