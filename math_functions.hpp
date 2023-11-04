#ifndef MATH_FUNCTIONS_HPP
#define MATH_FUNCTIONS_HPP

#include <vector>
#include <cmath>

double dejong1_function(const std::vector<double>& vec);
double schwefels_function(const std::vector<double>& vec);
double rastrigins_function(const std::vector<double>& vec);
double michalewiczs_function(const std::vector<double>& vec);
double h2p_function(const std::vector<double>& vec);

#endif