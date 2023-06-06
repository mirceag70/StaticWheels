#include "Helper.h"

void SoA_LP_gen_root_primes(void);

extern uint16_t firstPrimes[];
extern uint8_t root_primes_gap[];

uint8_t w8[] = { 1, 7, 11, 13, 17, 19, 23, 29 };

//w8_pattern[i][j] = (w8[i]*w8[j]) mod 30
uint8_t w8_pattern[8][8] = {{ 1,  7, 11, 13, 17, 19, 23, 29},
                            { 7, 19, 17,  1, 29, 13, 11, 23},
                            {11, 17,  1, 23,  7, 29, 13, 19},
                            {13,  1, 23, 19, 11,  7, 29, 17},
                            {17, 29,  7, 11, 19, 23,  1, 13},
                            {19, 13, 29,  7, 23,  1, 17, 11},
                            {23, 11, 13, 29,  1, 17, 19,  7},
                            {29, 23, 19, 17, 13, 11,  7,  1} };

//(w8[i] * w8_factors[i][j]) mod 30 = w8[j];
uint8_t w8_factors[8][8] = {{ 1,  7, 11, 13, 17, 19, 23, 29},
                            {13,  1, 23, 19, 11,  7, 29, 17},
                            {11, 17,  1, 23,  7, 29, 13, 19},
                            { 7, 19, 17,  1, 29, 13, 11, 23},
                            {23, 11, 13, 29,  1, 17, 19,  7},
                            {19, 13, 29,  7, 23,  1, 17, 11},
                            {17, 29,  7, 11, 19, 23,  1, 13},
                            {29, 23, 19, 17, 13, 11,  7,  1} };

constexpr unsigned sv8size = 8;
constexpr unsigned sv8max = 30;

uint8_t rsvR[NO_ROOT_PRIMES];

uint8_t vPrimes[1000];

void GetRootPrimes(void)
{
    SoA_LP_gen_root_primes();
    unsigned last_p = 5;
    for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
    {
        unsigned idx{};
        unsigned p = last_p + root_primes_gap[i];
        switch (p % 30)
        {
        case  1: idx = 0; break; case  7: idx = 1; break;
        case 11: idx = 2; break; case 13: idx = 3; break;
        case 17: idx = 4; break; case 19: idx = 5; break;
        case 23: idx = 6; break; case 29: idx = 7; break;
        default: __debugbreak();
        }
        rsvR[i] = idx;
        last_p = p;
    }
}

void SetBitInSv(uint64_t bitidx, uint8_t sv[])
{
    uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
    sv[idx] |= BIT_MASK[bit];
}

void SieveR(unsigned ridx, const uint64_t Nmax, const uint64_t svlen, uint8_t sv[])
{
    memset(sv, 0, svlen);
    uint64_t last_p = 5;
    for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
    {
        uint64_t p = last_p + root_primes_gap[i];
        if (p * p > Nmax) break;

        unsigned deltak = w8_factors[rsvR[i]][ridx];
        tpPrime k = deltak + 30ull * (p / 30);
        if (k < p) k += 30;
        for (; k * p <= Nmax; k += 30)
        {
            SetBitInSv(k * p / 30, sv);
        }

        last_p = p;
    }
}

uint64_t CountAllOnce(uint8_t vPrimes[], const uint64_t Nmax, const uint64_t svlen)
{
    cTimer tmr;
    tmr.Start();

    uint64_t numPrimes = 3; vPrimes[0] |= 1;
    for (unsigned i = 0; i < numPrimes; i++) AddPrime(firstPrimes[i]);
    uint64_t i;
    for (i = 0; i < (svlen-1); i++)
    {
        for (unsigned bit = 0; bit < 8; bit++)
        {
            uint64_t pbase = 30 * (8 * i + bit);
            for (unsigned s = 0; s < sv8size; s++)
            {
                uint8_t* sv = vPrimes + s * svlen;
                uint8_t flag = sv[i];
                if (!(flag & BIT_MASK[bit]))
                {
                    uint64_t p = pbase + w8[s];
                    numPrimes++;  AddPrime(p);
                }
            }
        }
    }
    //for (; i < svlen; i++)
    {
        for (unsigned bit = 0; bit < 8; bit++)
        {
            uint64_t pbase = 30 * (8 * i + bit);
            for (unsigned s = 0; s < sv8size; s++)
            {
                uint8_t* sv = vPrimes + s * svlen;
                uint8_t flag = sv[i];
                if (!(flag & BIT_MASK[bit]))
                {
                    uint64_t p = pbase + w8[s];
                    if (p > Nmax) break;
                    numPrimes++;  AddPrime(p);
                }
            }
        }
    }
    tmr.Stop(true, "counting");

    return numPrimes;
}

void CountOneChunk(uint8_t vPrimes[], const unsigned nt, const unsigned i, const uint64_t svlen, tpPrime* np)
{
    const tpPrime chunksz = svlen / nt;
    const tpPrime start = i * chunksz;
    const tpPrime stop = (i == (nt - 1)) ? (svlen - 1) : start + chunksz;
    *np = 0;

    for (uint64_t i = start; i < stop; i++)
    {
        for (unsigned bit = 0; bit < 8; bit++)
        {
            uint64_t pbase = 30 * (8 * i + bit);
            for (unsigned s = 0; s < sv8size; s++)
            {
                uint8_t* sv = vPrimes + s * svlen;
                uint8_t flag = sv[i];
                if (!(flag & BIT_MASK[bit]))
                {
                    uint64_t p = pbase + w8[s];
                    (*np)++;  AddPrimeThreaded(p);
                }
            }
        }
    }
}

uint64_t CountAllSegments(uint8_t vPrimes[], const uint64_t Nmax, const uint64_t svlen)
{
    cTimer tmr;
    tmr.Start();

    uint64_t numPrimes = 3; vPrimes[0] |= 1;
    for (unsigned i = 0; i < numPrimes; i++) AddPrime(firstPrimes[i]);
    
    const unsigned nt = 300; tpPrime np[nt]{ 0 };

    //count segments
    std::vector<std::thread> tc;
    for (auto i : Range<nt>())
        tc.push_back(std::thread(&CountOneChunk, vPrimes, nt, i, svlen, np + i));
    for (auto& t : tc) t.join();
    for (auto i : Range<nt>()) numPrimes += np[i];

    //count last bit
    for (unsigned bit = 0; bit < 8; bit++)
    {
        uint64_t pbase = 30 * (8 * (svlen-1) + bit);
        for (unsigned s = 0; s < sv8size; s++)
        {
            uint8_t* sv = vPrimes + s * svlen;
            uint8_t flag = sv[svlen - 1];
            if (!(flag & BIT_MASK[bit]))
            {
                uint64_t p = pbase + w8[s];
                if (p > Nmax) break;
                numPrimes++;  AddPrime(p);
            }
        }
    }

    tmr.Stop(true, "counting");

    return numPrimes;
}

uint64_t SoE_w8pat(const uint64_t Nmax, uint8_t vPrimes[], void*, void*)
{
    GetRootPrimes();

    const uint64_t Nsqrt = (uint64_t)floor(sqrt(Nmax));
    const uint64_t i_max = Nmax / (30ull * 8);
    const uint64_t svlen = i_max + 1;

    std::vector<std::thread> vt;
    for (unsigned i = 0; i < sv8size; i++)
    {
        vt.push_back(std::thread(&SieveR, i, Nmax, svlen, vPrimes + i * svlen));
        //SieveR(i, Nmax, svlen, vPrimes + i * svlen);
    }
    for (auto& t : vt) t.join();

    return CountAllSegments(vPrimes, Nmax, svlen);
    return CountAllOnce(vPrimes, Nmax, svlen);
}
