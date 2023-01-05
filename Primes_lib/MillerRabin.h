#include <array>
#include <cstdint>
#include <random>
#include <bitset>
#include <tuple>

#include "uint128_t-master/uint128_t.h"

#ifndef PRIMES_MILLER_RABIN_H
#define PRIMES_MILLER_RABIN_H

class MillerRabin {
private:
    constexpr const static std::array<int, 6> KNOWN_PRIMES = {11, 13, 17, 19, 23, 29};

    inline static std::random_device rd;
    inline static std::mt19937 mt19937 = std::mt19937(rd());
public:
    static uint32_t boundedRand(std::mt19937& rng, uint32_t range);
    static uint64_t boundedRand(uint64_t min, uint64_t max);

    static uint64_t ipow(uint64_t base, uint64_t exp);
    static uint128_t fastExp(uint128_t b, uint128_t e, uint128_t m, uint64_t k = 5);

    static long double li(uint64_t _x);
    static uint64_t pnt(uint64_t x);

    static bool trialComposite(
            uint64_t a,
            uint64_t s,
            uint64_t d,
            uint64_t n);

    static bool isPrime(uint64_t n, int t = 8);
};


#endif //PRIMES_MILLER_RABIN_H
