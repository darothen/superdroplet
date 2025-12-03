cmake_minimum_required(VERSION 3.20)
project(sd_cpp VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(
        -Wall
        -Wextra
        -Wpedantic
        -O3
        -march=native
        -ffast-math
        -flto
        -funroll-loops
        -fomit-frame-pointer
    )
    # Enable LTO for linking
    add_link_options(-flto)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_compile_options(/W4 /O2 /arch:AVX2)
endif()

# Source files
set(CORE_SOURCES
    src/core/droplet.cpp
)

set(PHYSICS_SOURCES
    src/physics/kernels.cpp
    src/physics/collision.cpp
)

set(UTILS_SOURCES
    src/utils/math.cpp
    src/utils/time.cpp
)

set(IO_SOURCES
    src/io/binning.cpp
)

# Create library
add_library(sd_cpp_lib
    ${CORE_SOURCES}
    ${PHYSICS_SOURCES}
    ${UTILS_SOURCES}
    ${IO_SOURCES}
)

target_include_directories(sd_cpp_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Create executable
add_executable(sd_cpp src/main.cpp)
target_link_libraries(sd_cpp PRIVATE sd_cpp_lib)

# Enable testing
enable_testing()
add_subdirectory(tests)

