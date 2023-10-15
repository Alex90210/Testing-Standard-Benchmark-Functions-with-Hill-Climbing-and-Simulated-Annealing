#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>

// hill climbing algorithm steps:
// 1. generate random bitstring
// 2. decode bitstring and evaluate it
// 3. do-while improve strategy (the best improvement, the least improvement, first improvement) on generated neighbourhood

unsigned L_binary_length(double interval_start, double interval_end, double epsilon) {
    unsigned pow_epsilon = 1 / epsilon;
    unsigned dim_number_of_bits = std::ceil(std::log2((interval_end - interval_start) * pow_epsilon));
    return dim_number_of_bits;
}

std::string generate_binary_string(double interval_start, double interval_end, double epsilon, unsigned number_of_dimensions) {
    unsigned pow_epsilon = 1 / epsilon;
    unsigned dim_number_of_bits = std::ceil(std::log2((interval_end - interval_start) * pow_epsilon));
    std::cout << dim_number_of_bits << std::endl;
    unsigned number_of_bits = number_of_dimensions * dim_number_of_bits;
    srand(time(0));
    std::string generated_string;
    for(size_t i {0}; i < number_of_bits; ++i) {
        generated_string += (rand() % 2 == 0) ? '1' : '0';
    }
    return generated_string;
}

std::vector<double> decode_binary_string(double interval_start, double interval_end, double epsilon, unsigned number_of_dimensions, std::string& binary_string) {
    // x will be between 0 and 2^n - 1, n is the length of the binary string
    std::vector <double> dimensional_values;
    unsigned dim_length = L_binary_length(interval_start, interval_end, epsilon);
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
    std::string fliped_string;
    for(auto i: binary_string) {
        if (i == '1')
            fliped_string += '0';
        else
            fliped_string += '1';
    }
    return fliped_string;
}

double de_jong_1(const std::vector<double>& vec) {
    double sum {0.0};
    for (auto x : vec)
        sum += x * x;
    return sum;
}

std::vector<double> generate_neighbourhood(double interval_start, double interval_end, double epsilon, unsigned number_of_dimensions, std::string& binary_string) {
    std::string flipped_binary_string = flip_bits(binary_string);
    std::vector<double> dim_values = decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, binary_string);
    // the de jung 1 value is needed for calculating the neighbours
    double function_value = de_jong_1(dim_values); // affects the modularity of the function, this needs to be changed
    unsigned interval_partitioner {2};
    std::vector<std::vector<double>> neighbours_in_all_dimensions;
    
}

int main () {

    double interval_start = -5.12;
    double interval_end = 5.12;
    double epsilon = 0.01;
    unsigned number_of_dimensions = 2;
    std::string test = generate_binary_string(interval_start, interval_end, epsilon, number_of_dimensions);
    std::cout << test << std::endl;
    std::vector<double> decoded = decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, test);
    for (auto i : decoded)
        std::cout << i << " ";
    std::cout << std::endl;
    std::cout << de_jong_1(decoded);
    /*std::string test = generate_binary_string(interval_start, interval_end, epsilon, number_of_dimensions);
    std::vector vec = decode_binary_string(interval_start, interval_end, epsilon, test, number_of_dimensions);
    std::cout << de_jong_1(vec) << std::endl;*/

    return 0;
}