#include "MillerRabin.h"
#include "uint128_t-master/uint128_t.h"

uint32_t MillerRabin::boundedRand(std::mt19937& rng, uint32_t range) {
    uint32_t x = rng();
    uint64_t m = uint64_t(x) * uint64_t(range);
    auto l = uint32_t(m);

    if (l < range) {
        uint32_t t = -range;
        if (t >= range) {
            t -= range;
            if (t >= range)
                t %= range;
        }
        while (l < t) {
            x = rng();
            m = uint64_t(x) * uint64_t(range);
            l = uint32_t(m);
        }
    }

    return m >> 32;
}

uint64_t MillerRabin::boundedRand(const uint64_t min, const uint64_t max) {
    return boundedRand(mt19937, max - min) + min;
}

uint64_t MillerRabin::ipow(uint64_t base, uint64_t exp) {
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
        case 255: // we use 255 as an overflow marker and return 0 on overflow/underflow
            if (base == 1) {
                return 1;
            }

            if (base == -1) {
                return 1 - 2 * (exp & 1);
            }

            return 0;
        //Do not lint duplicate branches, it is desired behavior to fallthrough through cases
        //NOLINTNEXTLINE(bugprone-branch-clone)
        case 6:
            if (exp & 1) result *= base;
            exp >>= 1;
            base *= base;
        case 5:
            if (exp & 1) result *= base;
            exp >>= 1;
            base *= base;
        case 4:
            if (exp & 1) result *= base;
            exp >>= 1;
            base *= base;
        case 3:
            if (exp & 1) result *= base;
            exp >>= 1;
            base *= base;
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

std::vector<uint128_t> precompute(uint128_t b, uint128_t m, uint64_t k) {
    uint128_t a = (b * b) % m;
    uint128_t s = 1 << k;

    std::vector<uint128_t> precomputed;
    precomputed.reserve(s);
    precomputed[0] = b % m;

    for (uint128_t i = 1; i < s; i++)
        precomputed[i] = (precomputed[i - 1] * a) % m;

    return precomputed;
}

std::tuple<uint32_t, uint32_t> splitD(uint128_t d) {
    uint32_t t = 0;
    while ((d & 0x1) == 0) {
        d >>= 1;
        t++;
    }

    return {t, d};
}

std::tuple<uint128_t, uint128_t> setupMask(uint128_t e, uint128_t k) {
    uint128_t mask = (1 << k) - 1;
    uint32_t p = ((32 + k-1)/k - 1) * k;
    mask <<= p;

    while (!(e & mask)) {
        mask >>= k;
        p -= k;
    }

    return {mask, p};
}

uint128_t MillerRabin::fastExp(uint128_t b, uint128_t e, uint128_t m, uint64_t k) {
    if (m == 1)
        return 0;
    if (b > m)
        b %= m;

    uint128_t r = 1;
    if (e < UINT32_MAX) {
        for (int64_t i = (uint32_t)1 << 29; i != 0; i >>= 1) {
            r = (r * r) % m;
            if (e & i)
                r = (r * b) % m;
        }
    } else {
        uint128_t B = 1 << k, d{}, t{}, o{};
        uint128_t x = 0;

        std::vector<uint128_t> precomputed = precompute(b, m, k);
        auto [mask, p] = setupMask(e, k);

        while (mask != 0) {
            d = (e & mask) >> p;
            x = (x * B) + d;

            if (d > 0) {
                const std::tuple<int32_t, int32_t> split = splitD(d);
                t = std::get<0>(split);
                o = std::get<1>(split);

                for (uint32_t i = 0; i < k-t; i++)
                    r = (r * r) % m;

                r = (r * precomputed[o/2]) % m;

                for (uint128_t i = 0; i < t; i++)
                    r = (r * r) % m;
            } else {
                for (uint128_t i = 0; i < k; i++)
                    r = (r * r) % m;
            }

            mask >>= k;
            p -= k;
        }
    }

    return r;
}

long double MillerRabin::li(const uint64_t _x) {
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

uint64_t MillerRabin::pnt(const uint64_t x) {
    return (uint64_t) li(x);
}

bool MillerRabin::trialComposite(
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

bool MillerRabin::isPrime(const uint64_t n, const int t) {
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