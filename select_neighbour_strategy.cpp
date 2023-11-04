#include "select_neighbour_strategy.hpp"

std::string first_improvement(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions,
                              const std::string& binary_string, double string_value,
                              double (*calculate_function)(const std::vector<double>& vec)) {

    double best_value {string_value};
    std::string copy_string = binary_string;


    for (unsigned i = get_random_unsigned(0, copy_string.length()); i < copy_string.length(); ++i) {

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

std::string best_improvement(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions,
                             const std::string& binary_string, double string_value, double (*calculate_function)(const std::vector<double>& vec)) {

    int index {-1};
    double best_value {string_value};
    std::string copy_string = binary_string;

    for (int i {0}; i < copy_string.length(); ++i) {

        copy_string[i] = (copy_string[i] == '1') ? '0' : '1'; // this could slow the process a lot, maybe a vector of bool values was a better idea
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

std::string best_improvement_up(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions,
                             const std::string& binary_string, double string_value, double (*calculate_function)(const std::vector<double>& vec)) {

    int index {-1};
    double best_value {string_value};
    std::string copy_string = binary_string;

    for (int i {0}; i < copy_string.length(); ++i) {

        copy_string[i] = (copy_string[i] == '1') ? '0' : '1'; // this could slow the process a lot, maybe a vector of bool values was a better idea
        double value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, copy_string));

        if (value > best_value) {
            best_value = value;
            index = i;
        }

        copy_string[i] = (copy_string[i] == '1') ? '0' : '1';
    }

    if (best_value > calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, copy_string)))
        copy_string[index] = (copy_string[index] == '1') ? '0' : '1';

    return copy_string;
}

bool descending_sort(double a, double b) {
    return a > b;
}

std::string worst_improvement(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions,
                              const std::string& binary_string, const double& string_value,
                              double (*calculate_function)(const std::vector<double>& vec)) {

    // I tried to solve it without storing the values, it did not work, original function is commented below
    std::string copy_string = binary_string;

    std::vector<double> value_vec;
    std::vector<std::string> string_vec;
    std::map<double, std::string> index_map;

    for (int i {0}; i < binary_string.length(); ++i) {

        copy_string[i] = (copy_string[i] == '1') ? '0' : '1';
        double value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, copy_string));

        value_vec.push_back(value);
        string_vec.push_back(copy_string);

        copy_string[i] = (copy_string[i] == '1') ? '0' : '1';
    }

    // mapping
    for (size_t i = 0; i < value_vec.size(); ++i) {
        index_map[value_vec[i]] = string_vec[i];
    }

    std::sort(value_vec.begin(), value_vec.end(), descending_sort);

    for (const auto& pair : index_map) {
        if(pair.first < string_value) {
            copy_string = pair.second;
            break;
        }
    }

    return copy_string;
}

// first attempt:
/*std::string best_improvement(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions,
                             const std::string& binary_string, double string_value, double (*calculate_function)(const std::vector<double>& vec)) {

    int index {-1};
    double best_value {string_value};
    double current_worst_improvement {-1000};
    std::string copy_string = binary_string;

    for (int i {0}; i < copy_string.length(); ++i) {

        copy_string[i] = (copy_string[i] == '1') ? '0' : '1';
        double value = calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, copy_string));

        if (value < best_value && value > current_worst_improvement) {
            current_worst_improvement = value;
            index = i;
        }

        copy_string[i] = (copy_string[i] == '1') ? '0' : '1';
    }

    if (current_worst_improvement < calculate_function(decode_binary_string(interval_start, interval_end, epsilon, number_of_dimensions, copy_string)))
        copy_string[index] = (copy_string[index] == '1') ? '0' : '1';

    return copy_string;
}*/

// second attempt:
/*std::string worst_improvement(const double& interval_start, const double& interval_end, double epsilon, unsigned number_of_dimensions,
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
}*/