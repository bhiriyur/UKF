cmake_minimum_required(VERSION 3.7)
project(UKF)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES src/ukf.cpp src/main.cpp)
add_executable(UKF ${SOURCE_FILES})