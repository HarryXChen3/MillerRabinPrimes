project(Primes)

set(HEADER_FILES
        MillerRabin.h
        uint128_t-master/uint128_t.h)

set(SOURCE_FILES
        MillerRabin.cpp
        uint128_t-master/uint128_t.cpp)

add_library(Primes_lib STATIC ${SOURCE_FILES} ${HEADER_FILES})