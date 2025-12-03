cmake_minimum_required(VERSION 3.20)

# Find GoogleTest package
find_package(GTest QUIET)

if(NOT GTest_FOUND)
    message(STATUS "GoogleTest not found, fetching from GitHub...")
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

# Enable testing
enable_testing()

# Helper macro to add tests
macro(add_sd_test test_name)
    add_executable(${test_name} ${test_name}.cpp)
    target_link_libraries(${test_name} PRIVATE sd_cpp_lib GTest::gtest GTest::gtest_main)
    target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    add_test(NAME ${test_name} COMMAND ${test_name})
endmacro()

# Add all test executables
add_sd_test(test_constants)
add_sd_test(test_droplet)
add_sd_test(test_math)
add_sd_test(test_time)
add_sd_test(test_kernels)
add_sd_test(test_collision)
add_sd_test(test_binning)

# Add a custom target to run all tests
add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure --verbose
    DEPENDS 
        test_constants 
        test_droplet 
        test_math 
        test_time 
        test_kernels 
        test_collision 
        test_binning
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
