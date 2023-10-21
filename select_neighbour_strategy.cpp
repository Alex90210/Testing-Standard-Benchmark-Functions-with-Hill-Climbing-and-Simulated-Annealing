#include "select_neighbour_strategy.hpp"

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