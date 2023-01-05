#include <algorithm>
#include <vector>

#include "gtest/gtest.h"
#include "MillerRabin.h"

TEST(MillerRabin, Ipow) {
    EXPECT_EQ(MillerRabin::ipow(1, 0), 1);
    EXPECT_EQ(MillerRabin::ipow(2, 6), 64);
    EXPECT_EQ(MillerRabin::ipow(3, 10), 59049);
    EXPECT_EQ(MillerRabin::ipow(6, 6), 46656);
    EXPECT_EQ(MillerRabin::ipow(2, 20), 1.048576E+6);
}

TEST(MillerRabin, FastExp) {
    EXPECT_EQ(MillerRabin::fastExp(2, 10, 145), 9);
    EXPECT_EQ(MillerRabin::fastExp(10, 780, 450), 100);
    EXPECT_EQ(MillerRabin::fastExp(6500, 8900, 999999000001), 558044759020);
    EXPECT_EQ(MillerRabin::fastExp(538570860, 15624984375, 999999000001), 787690165266);
}

TEST(MillerRabin, KnownPrimes) {
    EXPECT_EQ(MillerRabin::isPrime(7), true);
    EXPECT_EQ(MillerRabin::isPrime(8191), true);
    EXPECT_EQ(MillerRabin::isPrime(131071), true);
    EXPECT_EQ(MillerRabin::isPrime(524287), true);
    EXPECT_EQ(MillerRabin::isPrime(6700417), true);
    EXPECT_EQ(MillerRabin::isPrime(2147483647), true);
    //EXPECT_EQ(MillerRabin::isPrime(999999000001), true);
    //EXPECT_EQ(MillerRabin::isPrime(67280421310721), true);
}

TEST(MillerRabin, NPrimes) {
    const static uint128_t max = 10000;
    std::vector<uint128_t> primes;
    primes.reserve(MillerRabin::pnt(max));

    for (uint128_t i = 0; i < max; i++)
        if (MillerRabin::isPrime(i))
            primes.push_back(i);

    EXPECT_TRUE(std::binary_search(primes.begin(), primes.end(), 443));
    EXPECT_TRUE(std::binary_search(primes.begin(), primes.end(), 1901));
    EXPECT_TRUE(std::binary_search(primes.begin(), primes.end(), 7127));
    EXPECT_TRUE(std::binary_search(primes.begin(), primes.end(), 8191));
    EXPECT_TRUE(std::binary_search(primes.begin(), primes.end(), 9941));
}