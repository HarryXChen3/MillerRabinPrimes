#include <iostream>
#include <array>
#include <cstdint>
#include <random>
#include <chrono>

uint64_t boundedRand(std::mt19937_64& rng, uint64_t range) {
    uint64_t x = rng();
    if (range >= 1ull << 63) {
        while (x >= range)
            x = rng();
        return x;
    }

    uint64_t l = x * range;
    if (l < range) {
        uint64_t t = -range;
        t -= range;
        if (t >= range)
            t %= range;
        while (l < t) {
            x = rng();
            l = x * range;
        }
    }

    return l;
}

std::random_device rd;
std::mt19937_64 mt19937_64(rd());
uint64_t boundedRand(const uint64_t min, const uint64_t max) {
    return boundedRand(mt19937_64, max - min);
}

uint64_t ipow(uint64_t base, uint64_t exp) {
    static const std::array<uint8_t, 256> highest_bit_set = {
        0, 1, 2, 2, 3, 3, 3, 3,
        4, 4, 4, 4, 4, 4, 4, 4,
        5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5,
        6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 255, // anything past 63 is a guaranteed overflow with base > 1
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
    };

    uint64_t result = 1;
    switch (highest_bit_set[exp]) {
        case 255:
            if (base == 1)
                return 1;

            if (base == -1)
                return 1 - 2 * (exp & 1);

            return 0;
        case 6:
        case 5:
        case 4:
        case 3:
        case 2:
            if (exp & 1) result *= base;
            exp >>= 1;
            base *= base;
        case 1:
            if (exp & 1) result *= base;
        default:
            return result;
    }
}

uint64_t fastExp(uint64_t b, uint64_t e, uint64_t m) {
    if (m == 1)
        return 0;
    else {
        uint64_t r = 1;
        b %= m;
        while (e > 0) {
            if (e % 2 == 1)
                r = (r * b) % m;
            b = (b * b) % m;
            e >>= 1;
        }

        return r;
    }
}

bool trialComposite(
        const uint64_t a,
        const uint64_t d,
        const uint64_t s,
        const uint64_t n) {

    if (fastExp(a, d, n) == 1)
        return false;

    for (int i = 0; i < s; i++) {
        if (fastExp(a, d * ipow(2, i), n) == n - 1)
            return false;
    }

    return true;
}

const std::array<int, 6> KNOWN_PRIMES = {11, 13, 17, 19, 23, 29};
bool isPrime(const uint64_t n, const int t = 8) {
    if (n < 2 || n == 4 || n == 6 || n == 8 || n == 9)
        return false;

    if (n == 2 || n == 3 || n == 5 || n == 7)
        return true;

    for (const int prime : KNOWN_PRIMES)
        if (n % prime == 0)
            return n == prime;

    uint64_t s = 0;
    uint64_t d = n - 1;
    while (d % 2 == 0) {
        d >>= 1;
        s += 1;
    }

    for (int _i = 0; _i < t; _i++) {
        const uint64_t a = boundedRand(2, n - 1);

        if (trialComposite(a, d, s, n))
            return false;
    }

    return true;
}


long double li(const uint64_t _x) {
    static long double f = INT64_MIN;
    long double x = logl(_x);
    if (std::abs(x - 10) >= 12) {
        auto j = (uint64_t) (5.0l + (20l / std::abs(x)));
        while (j != 0) {
            f = (1.0l / ((1.0l / x) - (1.0l / j))) + x;
            j -= 1;
        }
        return expl(x) / f;
    } else if(x == 0)
        return f;
    else {
        auto j = (uint64_t) (10.0l + (2.0l * std::abs(x)));
        auto ij = j + 1;
        f = 1.0l / (ij * ij);

        while (j != 0) {
            f = ((f * j * x) + 1) / (j * j);
            j -= 1;
        }

        return (f * x) + logl(1.781072418l * std::abs(x));
    }
}

uint64_t pnt(const uint64_t x) {
    return (uint64_t) li(x);
}

void display(const std::vector<uint64_t>& primes) {
    const uint64_t size = primes.size();
    if (size > 20) {
        const std::vector<std::vector<int>> minifiedIndexes = {{0, 1, 2, 3, 4}, {5, 4, 3, 2, 1}};
        for (const std::vector<int>& indexes : minifiedIndexes) {
            if (minifiedIndexes.front() == indexes)
                for (int i = 0; i < indexes.size(); i++)
                    printf((i == indexes.size() - 1) ? "%llu..." : "%llu, ", primes.at(indexes.at(i)));
            else
                for (int i = 0; i < indexes.size(); i++)
                    printf((i == indexes.size() - 1) ? "%llu]\n" : "%llu, ", primes.at(size - indexes.at(i)));
        }
    } else {
        for (int i = 0; i < size; i++) {
            printf((i == size - 1) ? "%llu]\n" : "%llu, ", primes.at(i));
        }
    }
}

int main() {
    const uint64_t max = 7000000;

    std::vector<uint64_t> primes;
    primes.reserve(pnt(max));

    auto start = std::chrono::steady_clock::now();

    for (uint64_t i = 0; i < max; i++)
        if (isPrime(i))
            primes.push_back(i);

    auto end = std::chrono::steady_clock::now();
    auto diff = std::chrono::duration<double, std::milli>(end - start).count();

    printf("done in %.6fms\n", diff);
    display(primes);
}