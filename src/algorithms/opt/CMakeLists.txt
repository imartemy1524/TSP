cmake_minimum_required(VERSION 3.29)
project(opt)
set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++23")

add_executable(opt
        lin_kernighan.cpp)
