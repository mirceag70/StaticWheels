#include "Helper.h"

uint64_t SoE_6k(const uint64_t Nmax, uint8_t vPrimes[], void*, void*)
{
    auto idx2no = [](uint64_t idx) { return 3 * idx + 5 - (idx & 1); };
    auto no2idx = [](uint64_t no) { return no / 3 - 1; };
    auto SetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        vPrimes[idx] |= BIT_MASK[bit];
    };
    auto GetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        return (vPrimes[idx] & BIT_MASK[bit]);
    };

    const uint64_t i_max = Nmax / (3ull * 8);
    const uint64_t Nsqrt = (uint64_t)floor(sqrt(Nmax));
    unsigned numPrimes = 2; AddPrime(2); AddPrime(3);

    uint64_t i, j;
    for (i = 0; i <= i_max; i++) vPrimes[i] = 0;

    uint8_t stepi = 2;
    for (i = 5; i <= Nsqrt; i += stepi, stepi = 6 - stepi)
        if (!GetBit(no2idx(i)))
        {
            numPrimes++; AddPrime(i);

            uint64_t stepj2 = 2 * i, stepj6 = 6 * i;
            uint64_t stepj = ((i - (i / 3) * 3) == 1) ? (stepj6 - stepj2) : stepj2;

            for (j = i * i; j <= Nmax; j += stepj, stepj = stepj6 - stepj)
                SetBit(no2idx(j));
        }

    for ( /*continue counting primes*/; i <= Nmax; i += stepi, stepi = 6 - stepi)
        if (!GetBit(no2idx(i)))
        {
            numPrimes++; AddPrime(i);
        }

    return numPrimes;
}

void InitializeRootSieve2(const unsigned char N);

extern const unsigned rootN;
extern uint16_t firstPrimes[];
extern unsigned char root_sieve[];
extern unsigned root_sieve_size;// , root_sieve_N;

unsigned nextWheelIdx(unsigned idx) { return ((idx < (root_sieve_size - 1)) ? idx + 1 : 0); }

uint64_t SoE_rooted(const uint64_t Nmax, uint8_t vPrimes[], void*, void*)
{
    InitializeRootSieve2(rootN);

    //auto idx2no = [](uint64_t idx) { return 3 * idx + 5 - (idx & 1); };
    //auto no2idx = [](uint64_t no) { return no / 3 - 1; };
    auto SetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        vPrimes[idx] |= BIT_MASK[bit];
    };
    auto GetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        return (vPrimes[idx] & BIT_MASK[bit]);
    };

    const uint64_t i_max = Nmax / (3ull * 8);
    const uint64_t Nsqrt = (uint64_t)floor(sqrt(Nmax));
    unsigned numPrimes = rootN;
    for (unsigned i = 0; i < rootN; i++) AddPrime(firstPrimes[i]);

    for (uint64_t i = 0; i <= i_max; i++) vPrimes[i] = 0;

    uint64_t p = 1ull + root_sieve[0];
    unsigned ip = 1;
    for (; p <= Nsqrt; p += root_sieve[ip], ip = nextWheelIdx(ip))
        if (!GetBit(no2idx(p)))
        {
            numPrimes++; AddPrime(p);
            uint64_t n = 1ull + root_sieve[0];
            unsigned in = 1;
            for (; n*p <= Nmax; n += root_sieve[in], in = nextWheelIdx(in))
                SetBit(no2idx(n*p));
        }

    for (/*continue counting primes*/; p <= Nmax; p += root_sieve[ip], ip = nextWheelIdx(ip))
        if (!GetBit(no2idx(p)))
        {
            numPrimes++; AddPrime(p);
        }

    return numPrimes;
}

uint64_t SoE_rooted1(const uint64_t Nmax, uint8_t vPrimes[], void*, void*)
{
    InitializeRootSieve2(rootN);

    //auto idx2no = [](uint64_t idx) { return 3 * idx + 5 - (idx & 1); };
    //auto no2idx = [](uint64_t no) { return no / 3 - 1; };
    auto SetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        vPrimes[idx] |= BIT_MASK[bit];
    };
    auto GetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        return (vPrimes[idx] & BIT_MASK[bit]);
    };

    const uint64_t i_max = Nmax / (3ull * 8);
    const uint64_t Nsqrt = (uint64_t)floor(sqrt(Nmax));
    unsigned numPrimes = rootN;
    for (unsigned i = 0; i < rootN; i++) AddPrime(firstPrimes[i]);

    for (uint64_t i = 0; i <= i_max; i++) vPrimes[i] = 0;

    uint64_t p = 1ull + root_sieve[0];
    unsigned ip = 1;
    for (; p <= Nsqrt; p += root_sieve[ip], ip = nextWheelIdx(ip))
        if (!GetBit(no2idx(p)))
        {
            numPrimes++; AddPrime(p);
            uint64_t n = 1ull + root_sieve[0];
            unsigned in = 1;
            for (; n * p <= Nmax; n += root_sieve[in], in = nextWheelIdx(in))
                SetBit(no2idx(n * p));
        }

    for (/*continue counting primes*/; p <= Nmax; p += root_sieve[ip], ip = nextWheelIdx(ip))
        if (!GetBit(no2idx(p)))
        {
            numPrimes++; AddPrime(p);
        }

    return numPrimes;
}
