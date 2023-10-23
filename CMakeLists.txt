cmake_minimum_required(VERSION 3.26)
project(Homework1)

set(CMAKE_CXX_STANDARD 17)

add_executable(Homework1 main.cpp
        math_functions.hpp
        math_functions.cpp
        binary_operations.hpp
        binary_operations.cpp
        select_neighbour_strategy.hpp
        select_neighbour_strategy.cpp
        algorithms.hpp
        algorithms.cpp)

