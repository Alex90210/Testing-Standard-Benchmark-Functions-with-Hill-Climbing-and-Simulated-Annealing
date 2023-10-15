#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>

// hill climbing algorithm steps:
// 1. generate random bitstring
// 2. decode bitstring and evaluate it
// 3. do-while improve strategy (the best improvement, the least improvement, first improvement) on generated neighbourhood

unsigned D_binary_length(double interval_start, double interval_end, double epsilon) {
    unsigned pow_epsilon = 1 / epsilon;
    unsigned dim_number_of_bits = std::ceil(std::log2((interval_end - interval_start) * pow_epsilon));
    return dim_number_of_bits;
}

std::string generate_binary_string(double interval_start, double interval_end, double epsilon, unsigned number_of_dimensions) {
    unsigned pow_epsilon = 1 / epsilon;
    unsigned dim_number_of_bits = std::ceil(std::log2((interval_end - interval_start) * pow_epsilon));
    // std::cout << dim_number_of_bits << std::endl;
    unsigned number_of_bits = number_of_dimensions * dim_number_of_bits;
    srand(time(0));
    std::string generated_string;
    for(size_t i {0}; i < number_of_bits; ++i) {
        generated_string += (rand() % 2 == 0) ? '1' : '0';
    }
    return generated_string;
}

std::vector<double> decode_binary_string(double interval_start, double interval_end, double epsilon, unsigned number_of_dimensions, const std::string& binary_string) {
    // x will be between 0 and 2^n - 1, n is the length of the binary string
    std::vector <double> dimensional_values;
    unsigned dim_length = D_binary_length(interval_start, interval_end, epsilon);
    for (size_t i {0}; i < number_of_dimensions; ++i) {
        int xb_value {0};
        for (size_t j = dim_length * i; j < dim_length * (i + 1); ++j) {
            xb_value *= 2;
            xb_value += binary_string[j] - '0';
        }
        double x_value = xb_value / (pow(2, dim_length) - 1);
        x_value *= (interval_end - interval_start);
        x_value += interval_start; // x is now between interval_start and interval_end
        dimensional_values.push_back(x_value);
    }
    return dimensional_values;
}

std::string flip_bits(const std::string& binary_string) {
    std::string flipped_string;
    for(auto i: binary_string) {
        if (i == '1')
            flipped_string += '0';
        else
            flipped_string += '1';
    }
    return flipped_string;
}

double de_jong_1(const std::vector<double>& vec) {
    double sum {0.0};
    for (auto x : vec)
        sum += x * x;
    return sum;
}

std::vector<double> generate_neighbourhood(double interval_start, double interval_end, double epsilon, unsigned number_of_dimensions, const std::string& binary_string) {
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
}

double best_improvement(const std::vector<double>& vec, const double& best_solution) {
    double solution {best_solution};
    for (auto i : vec)
        if(i < solution)
            solution = i;

    return solution;
}

double hill_climbing(double interval_start, double interval_end, double epsilon,
                   unsigned number_of_dimensions, unsigned iterations) {
    double best_solution = 100000000;
    bool is_local_minimum = false;
    // generating initial solution
    std::string binary_string = generate_binary_string(interval_start, interval_end, epsilon, number_of_dimensions);
    std::vector<double> initial_vec = decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, binary_string);
    double initial_solution = de_jong_1(initial_vec);

    bool local_minimum {false};
    for (size_t i {0}; i < iterations; ++i) {
        std::vector<double> neighbours = generate_neighbourhood(interval_start, interval_end, epsilon, number_of_dimensions, binary_string);
        double best_option = best_improvement(neighbours, initial_solution);
        if (initial_solution == best_option)
            break;
        else if (best_option < initial_solution && best_option < best_solution)
            best_solution = best_option;
        // ............ 2 loops needed, to be continued
    }




    return best_solution;
    /*do {
        bool is_local_minimum = false;
        // generating initial solution
        std::string binary_string = generate_binary_string(interval_start, interval_end, epsilon, number_of_dimensions);
        std::vector<double> initial_vec = decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, binary_string);
        double solution = de_jong_1(initial_vec);

        do {
            std::vector<double> neighbours = generate_neighbourhood(interval_start, interval_end, epsilon, number_of_dimensions, binary_string);
            double this_solution = best_improvement(neighbours, solution);
            if (this_solution < solution)
                solution = this_solution;
            else is_local_minimum = true;
        } while (!is_local_minimum);

        --iterations;
        if (solution < best_solution)
            best_solution = solution;

    } while (iterations > 0);*/
} // hill climbing probably flawed, my best solution would be in the neighborhood of a random location

int main () {

    double interval_start = -5.12;
    double interval_end = 5.12;
    double epsilon = 0.01;
    unsigned number_of_dimensions = 2;
    unsigned iterations {100};




    double best = hill_climbing(interval_start, interval_end, epsilon, number_of_dimensions, iterations);
    std::cout << best;









    /*std::string test = generate_binary_string(interval_start, interval_end, epsilon, number_of_dimensions);

    std::vector<double> vec = decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, test);
    for (auto i : vec)
        std::cout << i << " ";
    std::cout << std::endl << "Valoarea initiala: " << de_jong_1(vec) << std::endl;

    generate_neighbourhood(interval_start, interval_end, epsilon, number_of_dimensions, test);
    std::vector<double> neighbourhood = generate_neighbourhood(interval_start, interval_end, epsilon, number_of_dimensions, test);

    for(auto i: neighbourhood)
        std::cout << i << std::endl;*/

    return 0;
}