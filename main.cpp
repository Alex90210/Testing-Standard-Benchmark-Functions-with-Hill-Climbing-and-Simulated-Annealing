#include "algorithms.hpp"
#include <chrono>
#include <iomanip>

int main () {

    std::string mode {"BI"};
    unsigned number_of_dimensions {30};
    double epsilon {0.001};
    unsigned iterations {100};
    double interval_start {-5.12};
    double interval_end {5.12};

    double temperature {100};

    auto start = std::chrono::high_resolution_clock::now();


    interval_start = -500;
    interval_end = 500;
    double sim = simulated_annealing_with_linear_decay(interval_start, interval_end, epsilon, number_of_dimensions, iterations, temperature, schwefels_function);
    std::cout << "Simulated annealer: " << sim << std::endl;
    double hill = hill_climbing(interval_start, interval_end, epsilon, number_of_dimensions, iterations, mode, schwefels_function);
    std::cout << "Hill climber: " << hill << std::endl;

    double random_point = rastrigins_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, generate_binary_string(interval_start, interval_end, epsilon, number_of_dimensions)));
    std::cout << "random rastrigin value: " << random_point << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Program executed in: " << duration.count() << " seconds." << std::endl;

    return 0;
}


  // main code
     // Schwefel
    // must be under -10000 (for 30 dimensions)
/*
    interval_start = -500;
    interval_end = 500;
    double best_s = hill_climbing(interval_start, interval_end, epsilon, number_of_dimensions, iterations, mode, schwefels_function);
    std::cout << std::fixed << std::setprecision(5) << "Schwefel: " << best_s << std::endl;
*/
    // Rastrigin
    // must be under -25 (for 30 dimensions)

    /*interval_start = -5.12;
    interval_end = 5.12;
    float best_r = hill_climbing(interval_start, interval_end, epsilon, number_of_dimensions, iterations, mode, rastrigins_function);
    std::cout << std::fixed << std::setprecision(5) << "Rastrigin: " << best_r << std::endl;*/

    /*// Michaleiwcz
    // must be under -25 (for 30 dimensions)

    interval_start = 0;
    interval_end = M_PI;
    double best_m = hill_climbing(interval_start, interval_end, epsilon, number_of_dimensions, iterations, mode, michalewiczs_function);
    std::cout << std::fixed << std::setprecision(5) << "Michalewicz: " << best_m << std::endl;*/

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