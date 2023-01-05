cmake_minimum_required(VERSION 3.23)
project(Primes)

set(CMAKE_CXX_STANDARD 20)

set(SOURCE_FILES main.cpp)
add_executable(Primes_run ${SOURCE_FILES})

include_directories(Primes_lib)
add_subdirectory(Primes_lib)

target_link_libraries(Primes_run Primes_lib)

add_subdirectory(tests)