cmake_minimum_required(VERSION 3.3)
project(sd_cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    droplet.cpp
    droplet.hpp
    sce.cpp
    util.cpp
    util.hpp constants.hpp)

add_executable(sd_cpp ${SOURCE_FILES})