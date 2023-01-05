project(tests)
add_subdirectory(googletest-release-1.12.1)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(Google_Tests_run KnownPrimesTest.cpp)
target_link_libraries(Google_Tests_run Primes_lib)
target_link_libraries(Google_Tests_run gtest gtest_main)