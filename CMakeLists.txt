cmake_minimum_required(VERSION 3.16)
project(numfft CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_library(numfft src/fft.cpp)
target_include_directories(numfft PUBLIC include)

enable_testing()
add_executable(test_fft tests/test_fft.cpp)
target_link_libraries(test_fft PRIVATE numfft)
add_test(NAME fft_tests COMMAND test_fft)

add_executable(bench_fft bench/bench.cpp)
target_link_libraries(bench_fft PRIVATE numfft)
