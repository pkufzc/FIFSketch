cmake_minimum_required(VERSION 3.9)
project(Bench)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-mcmodel=large -O2")

add_executable(bench main.cpp)
