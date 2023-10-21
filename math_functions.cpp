#include "math_functions.hpp"

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