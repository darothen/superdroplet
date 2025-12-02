cmake_minimum_required(VERSION 3.1)
project(sd_cpp)

#set( GCC-VERSION "4.8")
#if( GCC-VERSION )
#    set(CMAKE_C_COMPILER "gcc-${GCC-VERSION}")
#    set(CMAKE_CXX_COMPILER "g++-${GCC-VERSION}")
#endif( GCC-VERSION )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

set(SOURCE_FILES
    droplet.cpp
    droplet.hpp
    sce.cpp
    util.cpp
    util.hpp collisions.cpp collisions.hpp)

find_package( Boost COMPONENTS format random )
include_directories( ${Boost_INCLUDE_DIRS} )

find_package(OpenMP)

add_executable(sce ${SOURCE_FILES})
target_link_libraries(sce PUBLIC OpenMP::OpenMP_CXX)
