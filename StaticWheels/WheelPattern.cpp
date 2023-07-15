#include "Helper.h"

void SoA_LP_gen_root_primes(void);

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

uint8_t vPrimesx[1000];

void GetRootPrimes(void)
{
    SoA_LP_gen_root_primes();
    unsigned last_p = root_primes_gap[0];
    for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
    {
        unsigned idx{};
#if LIMIT_lg10 < 9
        unsigned p = last_p + root_primes_gap[i];
#else
        unsigned p = last_p + root_primes_gap[i] * 2;
#endif
        switch (p % 30)
        {
        default: NEVERHERE;
        case  1: idx = 0; break; case  7: idx = 1; break;
        case 11: idx = 2; break; case 13: idx = 3; break;
        case 17: idx = 4; break; case 19: idx = 5; break;
        case 23: idx = 6; break; case 29: idx = 7; break;
        }
        rsvR[i] = idx;
        last_p = p;
    }
}

void SieveR(unsigned ridx, const uint64_t stop, const uint64_t svlen, uint8_t sv[])
{
    memset(sv, 0, svlen);
    uint64_t last_p = 5;
    for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
    {
        uint64_t p = last_p + root_primes_gap[i];
        if (p * p > stop) break;

        unsigned deltak = w8_factors[rsvR[i]][ridx];
        tpPrime k = deltak + 30ull * (p / 30);
        if (k < p) k += 30;
        for (; k * p <= stop; k += 30)
        {
            SetBit(k * p / 30, sv);
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


uint64_t CSA(uint64_t& a, uint64_t  b, uint64_t c)
{
    uint64_t v = a & b | (a ^ b) & c;
    a ^= b ^ c;
    return v;
}
uint64_t OnesCount(uint64_t x)
{
    x -= (x >> 1) & 0x5555'5555'5555'5555;
    x = ((x >> 2) & 0x3333'3333'3333'3333) + (x & 0x3333'3333'3333'3333);
    return (uint64_t)((((x + (x >> 4)) & 0x0F0F'0F0F'0F0F'0F0F) * 0x0101'0101'0101'0101) >> 56);
}
uint64_t OnesCount(uint8_t* p, uint64_t bits)
{
    uint64_t z, y, x;
    uint64_t c = 0;

    for (x = y = z = 0; bits >= 8; bits -= 8)
    {
        c += OnesCount(CSA(x, CSA(y, CSA(z, *p++, *p++),
                                     CSA(z, *p++, *p++)),
                              CSA(y, CSA(z, *p++, *p++),
                                     CSA(z, *p++, *p++))));
    }
    c = (c << 3) + (OnesCount(x) << 2) + (OnesCount(y) << 1) + OnesCount(z);

    for (; bits > 0; bits--)
        c += OnesCount(*p++);

    return c;
}
uint8_t OnesCount8(uint8_t x)
{
    x -= (x >> 1) & 0x55;
    x = ((x >> 2) & 0x33) + (x & 0x33);
    return (x >> 4) + (x & 0x0F);
}
constexpr uint64_t BLEN = 4'000'000'456;
void testCSA(void)
{
    uint8_t *bfr = new uint8_t[BLEN];

    for (unsigned i = 0; i < BLEN; i++) bfr[i] = i % 256;

    tpPrime bits1 = 0, bits2 = 0;
    cTimer tmr; tmr.Start();

    for (unsigned i = 0; i < BLEN; i++)
    {
        bits1 += bit0_count_table[bfr[i]];
        //bits2 += OnesCount8(bfr[i]);
        //assert(bit0_count_table[bfr[i]] == (8 - OnesCount8(bfr[i])));
    }
    bits2 = BLEN*8 - bits2;
    //bits2 = BLEN*8 - OnesCount(bfr, BLEN);

    tmr.Stop(true, "popcount");
    std::cout << bits1 << "|" << bits2;

    delete[] bfr;
}


void CountOneChunk(uint8_t vPrimes[], const unsigned nt, const unsigned i, const uint64_t svlen, tpPrime* np)
{
    auto do1bit = [&](uint8_t flg, uint8_t msk, uint64_t p)
    {
        if (!(flg & msk)) { AddPrimeThreaded(p); }
    };

    const tpPrime chunksz = svlen / nt;
    const tpPrime start = i * chunksz;
    const tpPrime stop = (i == (nt - 1)) ? (svlen - 1) : start + chunksz;
    tpPrime cnt = 0;

    uint8_t* pflg1 = vPrimes, * pflg7 = vPrimes + svlen, * pflg11 = vPrimes + 2 * svlen;
    uint8_t* pflg13 = vPrimes + 3 * svlen, * pflg17 = vPrimes + 4 * svlen, * pflg19 = vPrimes + 5 * svlen;
    uint8_t* pflg23 = vPrimes + 6 * svlen, * pflg29 = vPrimes + 7 * svlen;


    for (uint64_t i = start; i < stop; i++)
    {
        //uint64_t pbase;
        uint8_t flg1 = pflg1[i], flg7 = pflg7[i], flg11 = pflg11[i], flg13 = pflg13[i],
            flg17 = pflg17[i], flg19 = pflg19[i], flg23 = pflg23[i], flg29 = pflg29[i];
        cnt += (tpPrime)bit0_count_table[flg1] + bit0_count_table[flg7] + bit0_count_table[flg11] + bit0_count_table[flg13]
            + bit0_count_table[flg17] + bit0_count_table[flg19] + bit0_count_table[flg23] + bit0_count_table[flg29];
#ifdef GENERATE_PRIMES
        //effectively generate the primes in order
        {
            constexpr unsigned bit = 0; constexpr auto msk = BIT_MASK[bit];
            uint64_t pbase = 30 * (8 * i + bit);
            do1bit(flg1, msk, pbase + w8[0]); do1bit(flg7, msk, pbase + w8[1]);
            do1bit(flg11, msk, pbase + w8[2]); do1bit(flg13, msk, pbase + w8[3]);
            do1bit(flg17, msk, pbase + w8[4]); do1bit(flg19, msk, pbase + w8[5]);
            do1bit(flg23, msk, pbase + w8[6]); do1bit(flg29, msk, pbase + w8[7]);
            //if (!(flg1 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[0]); }
            //if (!(flg7 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[1]); }
            //if (!(flg11 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[2]); }
            //if (!(flg13 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[3]); }
            //if (!(flg17 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[4]); }
            //if (!(flg19 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[5]); }
            //if (!(flg23 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[6]); }
            //if (!(flg29 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[7]); }
        }
        {
            constexpr unsigned bit = 1; constexpr auto msk = BIT_MASK[bit];
            uint64_t pbase = 30 * (8 * i + bit);
            do1bit(flg1, msk, pbase + w8[0]); do1bit(flg7, msk, pbase + w8[1]);
            do1bit(flg11, msk, pbase + w8[2]); do1bit(flg13, msk, pbase + w8[3]);
            do1bit(flg17, msk, pbase + w8[4]); do1bit(flg19, msk, pbase + w8[5]);
            do1bit(flg23, msk, pbase + w8[6]); do1bit(flg29, msk, pbase + w8[7]);
        }
        {
            constexpr unsigned bit = 2; constexpr auto msk = BIT_MASK[bit];
            uint64_t pbase = 30 * (8 * i + bit);
            do1bit(flg1, msk, pbase + w8[0]); do1bit(flg7, msk, pbase + w8[1]);
            do1bit(flg11, msk, pbase + w8[2]); do1bit(flg13, msk, pbase + w8[3]);
            do1bit(flg17, msk, pbase + w8[4]); do1bit(flg19, msk, pbase + w8[5]);
            do1bit(flg23, msk, pbase + w8[6]); do1bit(flg29, msk, pbase + w8[7]);
        }
        {
            constexpr unsigned bit = 3; constexpr auto msk = BIT_MASK[bit];
            uint64_t pbase = 30 * (8 * i + bit);
            do1bit(flg1, msk, pbase + w8[0]); do1bit(flg7, msk, pbase + w8[1]);
            do1bit(flg11, msk, pbase + w8[2]); do1bit(flg13, msk, pbase + w8[3]);
            do1bit(flg17, msk, pbase + w8[4]); do1bit(flg19, msk, pbase + w8[5]);
            do1bit(flg23, msk, pbase + w8[6]); do1bit(flg29, msk, pbase + w8[7]);
        }
        {
            constexpr unsigned bit = 4; constexpr auto msk = BIT_MASK[bit];
            uint64_t pbase = 30 * (8 * i + bit);
            do1bit(flg1, msk, pbase + w8[0]); do1bit(flg7, msk, pbase + w8[1]);
            do1bit(flg11, msk, pbase + w8[2]); do1bit(flg13, msk, pbase + w8[3]);
            do1bit(flg17, msk, pbase + w8[4]); do1bit(flg19, msk, pbase + w8[5]);
            do1bit(flg23, msk, pbase + w8[6]); do1bit(flg29, msk, pbase + w8[7]);
        }
        {
            constexpr unsigned bit = 5; constexpr auto msk = BIT_MASK[bit];
            uint64_t pbase = 30 * (8 * i + bit);
            do1bit(flg1, msk, pbase + w8[0]); do1bit(flg7, msk, pbase + w8[1]);
            do1bit(flg11, msk, pbase + w8[2]); do1bit(flg13, msk, pbase + w8[3]);
            do1bit(flg17, msk, pbase + w8[4]); do1bit(flg19, msk, pbase + w8[5]);
            do1bit(flg23, msk, pbase + w8[6]); do1bit(flg29, msk, pbase + w8[7]);
        }
        {
            constexpr unsigned bit = 6; constexpr auto msk = BIT_MASK[bit];
            uint64_t pbase = 30 * (8 * i + bit);
            do1bit(flg1, msk, pbase + w8[0]); do1bit(flg7, msk, pbase + w8[1]);
            do1bit(flg11, msk, pbase + w8[2]); do1bit(flg13, msk, pbase + w8[3]);
            do1bit(flg17, msk, pbase + w8[4]); do1bit(flg19, msk, pbase + w8[5]);
            do1bit(flg23, msk, pbase + w8[6]); do1bit(flg29, msk, pbase + w8[7]);
        }
        {
            constexpr unsigned bit = 7; constexpr auto msk = BIT_MASK[bit];
            uint64_t pbase = 30 * (8 * i + bit);
            do1bit(flg1, msk, pbase + w8[0]); do1bit(flg7, msk, pbase + w8[1]);
            do1bit(flg11, msk, pbase + w8[2]); do1bit(flg13, msk, pbase + w8[3]);
            do1bit(flg17, msk, pbase + w8[4]); do1bit(flg19, msk, pbase + w8[5]);
            do1bit(flg23, msk, pbase + w8[6]); do1bit(flg29, msk, pbase + w8[7]);
        }
#endif
    }
    (*np) = cnt;
}

void CountOneChunkX(uint8_t vPrimes[], const unsigned nt, const unsigned i, const uint64_t svlen, tpPrime* np)
{
    const tpPrime chunksz = svlen / nt;
    const tpPrime start = i * chunksz;
    const tpPrime stop = (i == (nt - 1)) ? (svlen - 1) : start + chunksz;
    *np = 0;

    //uint8_t* pflg1 = vPrimes, * pflg7 = vPrimes + svlen, * pflg11 = vPrimes + 2 * svlen;
    //uint8_t* pflg13 = vPrimes + 3 * svlen, * pflg17 = vPrimes + 4 * svlen, * pflg19 = vPrimes + 5 * svlen;
    //uint8_t* pflg23 = vPrimes + 6 * svlen, * pflg29 = vPrimes + 7 * svlen;

    //for (uint64_t i = start; i < stop; i++)
    //{
    //    //uint64_t pbase = 30 * 8 * i;
    //    uint8_t flg1 = pflg1[i], flg7 = pflg7[i], flg11 = pflg11[i], flg13 = pflg13[i],
    //    flg17 = pflg17[i], flg19 = pflg19[i], flg23 = pflg23[i], flg29 = pflg29[i];
    //    {
    //        constexpr unsigned bit = 0;
    //        uint64_t pbase = 30 * (8 * i + bit);
    //        if (!(flg1 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[0]); }
    //        if (!(flg7 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[1]); }
    //        if (!(flg11 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[2]); }
    //        if (!(flg13 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[3]); }
    //        if (!(flg17 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[4]); }
    //        if (!(flg19 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[5]); }
    //        if (!(flg23 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[6]); }
    //        if (!(flg29 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[7]); }
    //    }
    //    {
    //        constexpr unsigned bit = 1;
    //        uint64_t pbase = 30 * (8 * i + bit);
    //        //pbase += 30;
    //        if (!(flg1 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[0]); }
    //        if (!(flg7 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[1]); }
    //        if (!(flg11 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[2]); }
    //        if (!(flg13 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[3]); }
    //        if (!(flg17 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[4]); }
    //        if (!(flg19 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[5]); }
    //        if (!(flg23 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[6]); }
    //        if (!(flg29 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[7]); }
    //    }
    //    {
    //        constexpr unsigned bit = 2;
    //        uint64_t pbase = 30 * (8 * i + bit);
    //        //pbase += 30;
    //        if (!(flg1 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[0]); }
    //        if (!(flg7 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[1]); }
    //        if (!(flg11 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[2]); }
    //        if (!(flg13 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[3]); }
    //        if (!(flg17 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[4]); }
    //        if (!(flg19 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[5]); }
    //        if (!(flg23 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[6]); }
    //        if (!(flg29 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[7]); }
    //    }
    //    {
    //        constexpr unsigned bit = 3;
    //        uint64_t pbase = 30 * (8 * i + bit);
    //        //pbase += 30;
    //        if (!(flg1 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[0]); }
    //        if (!(flg7 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[1]); }
    //        if (!(flg11 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[2]); }
    //        if (!(flg13 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[3]); }
    //        if (!(flg17 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[4]); }
    //        if (!(flg19 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[5]); }
    //        if (!(flg23 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[6]); }
    //        if (!(flg29 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[7]); }
    //    }
    //    {
    //        constexpr unsigned bit = 4;
    //        uint64_t pbase = 30 * (8 * i + bit);
    //        //pbase += 30;
    //        if (!(flg1 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[0]); }
    //        if (!(flg7 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[1]); }
    //        if (!(flg11 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[2]); }
    //        if (!(flg13 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[3]); }
    //        if (!(flg17 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[4]); }
    //        if (!(flg19 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[5]); }
    //        if (!(flg23 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[6]); }
    //        if (!(flg29 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[7]); }
    //    }
    //    {
    //        constexpr unsigned bit = 5;
    //        uint64_t pbase = 30 * (8 * i + bit);
    //        //pbase += 30;
    //        if (!(flg1 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[0]); }
    //        if (!(flg7 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[1]); }
    //        if (!(flg11 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[2]); }
    //        if (!(flg13 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[3]); }
    //        if (!(flg17 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[4]); }
    //        if (!(flg19 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[5]); }
    //        if (!(flg23 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[6]); }
    //        if (!(flg29 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[7]); }
    //    }
    //    {
    //        constexpr unsigned bit = 6;
    //        uint64_t pbase = 30 * (8 * i + bit);
    //        //pbase += 30;
    //        if (!(flg1 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[0]); }
    //        if (!(flg7 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[1]); }
    //        if (!(flg11 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[2]); }
    //        if (!(flg13 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[3]); }
    //        if (!(flg17 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[4]); }
    //        if (!(flg19 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[5]); }
    //        if (!(flg23 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[6]); }
    //        if (!(flg29 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[7]); }
    //    }
    //    {
    //        constexpr unsigned bit = 7;
    //        uint64_t pbase = 30 * (8 * i + bit);
    //        //pbase += 30;
    //        if (!(flg1 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[0]); }
    //        if (!(flg7 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[1]); }
    //        if (!(flg11 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[2]); }
    //        if (!(flg13 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[3]); }
    //        if (!(flg17 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[4]); }
    //        if (!(flg19 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[5]); }
    //        if (!(flg23 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[6]); }
    //        if (!(flg29 & BIT_MASK[bit])) { (*np)++;  AddPrimeThreaded(pbase + w8[7]); }
    //    }
    //}

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
    for (unsigned i = 0; i < numPrimes; i++) AddPrimeThreaded(firstPrimes[i]);
    
    const unsigned nt = 390; tpPrime np[nt]{ 0 };

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
                numPrimes++;  AddPrimeThreaded(p);
            }
        }
    }

    tmr.Stop(true, "counting");

    return numPrimes;
}

uint64_t SoE_w8pat(const uint64_t Nmax, uint8_t vPrimes[], void*, void*)
{
    GetRootPrimes();

    const uint64_t svlen = Nmax / (30ull * 8) + 1;

    std::vector<std::thread> vt;
    for (unsigned i = 0; i < sv8size; i++)
    {
        vt.push_back(std::thread(&SieveR, i, Nmax, svlen, vPrimes + i * svlen));
        //SieveR(i, Nmax, svlen, vPrimes + i * svlen);
    }
    for (auto& t : vt) t.join();

    return CountAllSegments(vPrimes, Nmax, svlen);
    //return CountAllOnce(vPrimes, Nmax, svlen);
}

//void SieveChunk()
//{
//
//}
//
//void SieveSegment(const tpPrime limit, uint8_t svall[], const unsigned ridx, const unsigned nt, const unsigned i)
//{
//    const uint64_t svlen = limit / (30ull * 8) + 1; // len for one r
//    uint8_t *svr = svall + ridx * svlen;            // sieve for current r
//    
//    const uint64_t chunksize = limit / (30ull * nt) ; // in bytes        
//    const uint64_t istart = i * chunksize, istop = istart + chunksize;
//    uint8_t* sv = svr + istart; memset(sv, 0, chunksize);
//    const uint64_t pstart = istart * 30, pstop = pstart * 30;
//
//    uint64_t last_p = 5;
//    for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
//    {
//        uint64_t p = last_p + root_primes_gap[i];
//        if (p * p > pstop) break;
//
//        unsigned deltak = w8_factors[rsvR[i]][ridx];
//        tpPrime k = deltak + 30ull * (p / 30);
//        if (k < p) k += 30;
//        for (; k * p <= pstop; k += 30)
//        {
//            SetBit(k * p / 30, svr);
//        }
//
//        last_p = p;
//    }
//}
//
//void ParallelSieve(const tpPrime limit, uint8_t sv[])
//{
//    const unsigned nt = 1;
//    const std::array<const unsigned, 8> vnt = 
//        { /*1*/nt + 0, /* 7*/nt + 0, /*11*/nt + 0, /*13*/nt + 0, 
//         /*17*/nt + 0, /*19*/nt + 0, /*23*/nt + 0, /*29*/nt + 0};
//    std::array<std::vector<std::thread>, 8> vt;
//
//    for (auto i : Range<8>())
//        for (unsigned j = 0; j < vnt[i]; j++)
//            vt[i].push_back(std::thread(&SieveSegment, limit, sv, i, vnt[i], j));
//
//    for (auto i : Range<8>()) for (auto& t : vt[i]) t.join();
//}
//
//uint64_t SoE_w8pat1(const uint64_t Nmax, uint8_t vPrimes[], void*, void*)
//{
//    GetRootPrimes();
//
//    ParallelSieve(Nmax, vPrimes);
//
//    const uint64_t svlen = Nmax / (30ull * 8) + 1;
//    return CountAllSegments(vPrimes, Nmax, svlen);
//    //return CountAllOnce(vPrimes, Nmax, svlen);
//}

constexpr tpPrime CHUNK_SIZE = 75'000'000ull;

void SieveOneChunkR(const unsigned ridx, const uint64_t start, const uint64_t stop, uint8_t sv[])
{
    //uint64_t last_p = 5, last_k = 30ull * (start / 7 / 30);
    //for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
    //{
    //    uint64_t p = last_p + root_primes_gap[i]; if (p * p > stop) break;
    //    unsigned deltak = w8_factors[rsvR[i]][ridx];
    //    tpPrime k = deltak + last_k; 
    //    tpPrime n = k * p, deltan = p * 30;
    //    if (k * p < start) 
    //        k += 30;
    //    else
    //    {
    //        while (k * p >= (start + p * 30)) k -= 30;
    //        if (k < p) k += 30; last_k = k - deltak;
    //    }
    //    tpPrime kk = deltak + 30ull * (start / p / 30); if (kk < p) kk += 30;
    //    if (kk * p < start) kk += 30;
    //    assert(k == kk);
    //    for (; k * p <= stop; k += 30)
    //    {
    //        SetBit(k * p / 30, sv);
    //    }
    //    last_p = p;
    //}
    uint64_t last_p = 5;
    for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
    {
        uint64_t p = last_p + root_primes_gap[i]; if (p * p > stop) break;
        unsigned deltak = w8_factors[rsvR[i]][ridx];
        tpPrime k = deltak + 30ull * (start / p / 30); if (k == 1) k = 31; 
        tpPrime n = k * p, deltan = 30 * p;
        if (n < start) n += deltan; assert(n >= start);
        for (; n <= stop; n += deltan)
        {
            SetBit(n / 30, sv);
        }
        last_p = p;
    }
}

void SieveChunksR(const unsigned ridx, const uint64_t stop, const uint64_t svlen, uint8_t sv[])
{
    memset(sv, 0, svlen);
    tpPrime iter_min = 0;
    if (stop > CHUNK_SIZE)
        for (; iter_min < (stop - CHUNK_SIZE); iter_min += CHUNK_SIZE)
        {
            SieveOneChunkR(ridx, iter_min, iter_min + CHUNK_SIZE, sv);
        }
    SieveOneChunkR(ridx, iter_min, stop, sv);
}

uint64_t SoE_w8pat_inc(const uint64_t Nmax, uint8_t vPrimes[], void*, void*)
{
    GetRootPrimes();

    const uint64_t svlen = Nmax / (30ull * 8) + 1;

    std::vector<std::thread> vt;
    for (unsigned i = 0; i < sv8size; i++)
    {
        vt.push_back(std::thread(&SieveChunksR, i, Nmax, svlen, vPrimes + i * svlen));
        //SieveChunksR(i, Nmax, svlen, vPrimes + i * svlen);
    }
    for (auto& t : vt) t.join();

    return CountAllSegments(vPrimes, Nmax, svlen);
    //return CountAllOnce(vPrimes, Nmax, svlen);
}

tpPrime lastns[8][NO_ROOT_PRIMES/10];

void SieveOneChunkRInc(const unsigned ridx, const uint64_t start, const uint64_t stop, uint8_t sv[])
{
    uint64_t last_p = 5;
    for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
    {
        const uint64_t p = last_p + root_primes_gap[i]; if (p * p > stop) break;
        const uint64_t deltan = 30 * p;
        //unsigned deltak = w8_factors[rsvR[i]][ridx];
        //tpPrime k = deltak + 30ull * (start / p / 30); if (k == 1) k = 31;
        //tpPrime n = k * p, deltan = 30 * p;
        tpPrime n = lastns[ridx][i];
        if (n < stop)
        {
            for (; n <= stop; n += deltan)
            {
                SetBit(n / 30, sv);
            }
            lastns[ridx][i] = n;
        }
        last_p = p;
    }
}

void SieveChunksRInc(const unsigned ridx, const uint64_t stop, const uint64_t svlen, uint8_t sv[])
{
    memset(sv, 0, svlen);
    //init multiples vector
    uint64_t p = 5;
    for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
    {
        p += root_primes_gap[i]; if (p * p > stop) break;
        const unsigned deltak = w8_factors[rsvR[i]][ridx];
        tpPrime k = deltak + 30ull * (p / 30); if (k < p) k += 30;
        const tpPrime n = k * p; lastns[ridx][i] = n;
    }

    tpPrime iter_min = 0;
    if (stop > CHUNK_SIZE)
        for (; iter_min < (stop - CHUNK_SIZE); iter_min += CHUNK_SIZE)
        {
            SieveOneChunkRInc(ridx, iter_min, iter_min + CHUNK_SIZE, sv);
        }
    SieveOneChunkRInc(ridx, iter_min, stop, sv);
}

uint64_t SoE_w8pat_inc1(const uint64_t Nmax, uint8_t vPrimes[], void*, void*)
{
    GetRootPrimes();

    const uint64_t svlen = Nmax / (30ull * 8) + 1;

    std::vector<std::thread> vt;
    for (unsigned i = 0; i < sv8size; i++)
    {
        //vt.push_back(std::thread(&SieveChunksRInc, i, Nmax, svlen, vPrimes + i * svlen));
        SieveChunksRInc(i, Nmax, svlen, vPrimes + i * svlen);
    }
    for (auto& t : vt) t.join();

    return CountAllSegments(vPrimes, Nmax, svlen);
    //return CountAllOnce(vPrimes, Nmax, svlen);
}

#undef CHUNK_SIZE
constexpr tpPrime CHUNK_SIZE_P = 75'000'000ull;

#undef lastns
thread_local tpPrime *lastnsp = NULL;

void SieveOneChunkRP(const unsigned ridx, const uint64_t start, const uint64_t stop, uint8_t sv[])
{
    const uint64_t maxn = stop / 30;
    const uint64_t maxp = (uint64_t)floor(sqrt(stop));
    uint64_t last_p = 5;
    for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
    {
        uint64_t p = last_p + root_primes_gap[i]; if (p > maxp) break;
        tpPrime n = lastnsp[i];
        if (n <= maxn)
        {
            //uint64_t idx0 = n / 8, bit0 = n % 8;
            //uint64_t idxdeltap = p / 8, idxdeltab = p % 8;
            for (; n <= maxn; n += p)
            {
                sv[n / 8] |= BIT_MASK[n % 8];

                //uint64_t idx = n / 8, bit = n % 8;
                //sv[idx] |= BIT_MASK[bit];
                //assert(idx0 == idx and bit0 == bit);
                //sv[idx0] |= BIT_MASK[bit0];
                //switch (idxdeltab)
                //{
                //default: NEVERHERE;
                //case 1: 
                //    if (bit0 == 7) { bit0 = 0; idx0++; } 
                //    else bit0++; break;
                //case 3: 
                //    if (bit0 < 5) { bit0 += 3; }
                //    else { bit0-=5; idx0++; }; break;
                //case 5:
                //    if (bit0 < 3) { bit0 += 5; }
                //    else { bit0 -= 3; idx0++; }; break;
                //case 7:
                //    if (bit0 == 0) { bit0 = 7; }
                //    else { bit0--; idx0++; }; break;
                //}
                //idx0 += idxdeltap;
            }
            lastnsp[i] = n;
        }
        last_p = p;
    }
}

void SieveChunksRP(const unsigned ridx, const uint64_t Nmax, uint8_t sv[], const unsigned nt, const unsigned i)
{
    //cTimer tmr; tmr.Start();

    const tpPrime chunksz = Nmax / nt;
    const tpPrime start = i * chunksz;
    const tpPrime stop = (i == (nt - 1)) ? Nmax: start + chunksz;

    //init multiples vector
    lastnsp = new tpPrime[NO_ROOT_PRIMES]; lastnsp[0] = 0;
    uint64_t p = 5;
    if (start == 0)
        for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
        {
            p += root_primes_gap[i]; if (p * p > stop) break;
            const unsigned deltak = w8_factors[rsvR[i]][ridx];
            tpPrime k = deltak;
                k += 30ull * (p / 30); if (k < p) k += 30;
            const tpPrime n = k * p; lastnsp[i] = n/30;
        }
    else
        for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
        {
            p += root_primes_gap[i]; if (p * p > stop) break;
            const unsigned deltak = w8_factors[rsvR[i]][ridx];
            tpPrime k = deltak;
                k += 30ull * (start / p / 30); if (k == 1) k = 31;
            const tpPrime n = k * p; lastnsp[i] = n / 30;
        }

    tpPrime iter_min = start;
    if (stop > start + CHUNK_SIZE_P)
        for (; iter_min < (stop - CHUNK_SIZE_P); iter_min += CHUNK_SIZE_P)
        {
            SieveOneChunkRP(ridx, iter_min, iter_min + CHUNK_SIZE_P, sv);
        }
    SieveOneChunkRP(ridx, iter_min, stop, sv);

    delete[] lastnsp;

    //std::cout << ridx << ":" << i << "|";
    //tmr.Stop(true, "sv");
}

void ParallelSieve(const tpPrime limit, uint8_t sv[], const uint64_t svlen)
{
    memset(sv, 0, 8*svlen);
    const unsigned nt = 3;
    const std::array<const unsigned, 8> vnt = 
        { /*1*/nt + 0, /* 7*/nt + 0, /*11*/nt + 0, /*13*/nt + 0, 
         /*17*/nt + 0, /*19*/nt + 0, /*23*/nt + 0, /*29*/nt + 0};
    std::array<std::vector<std::thread>, 8> vt;

    for (auto i : Range<8>())
        for (unsigned j = 0; j < vnt[i]; j++)
        {
            vt[i].push_back(std::thread(&SieveChunksRP, i, limit, sv + i * svlen, vnt[i], j));
            //SieveChunksRP( i, limit, sv + i * svlen, vnt[i], j);
        }
    for (auto i : Range<8>()) for (auto& t : vt[i]) t.join();
}

uint64_t SoE_w8pat_inc2(const uint64_t Nmax, uint8_t vPrimes[], void*, void*)
{
    GetRootPrimes();

    const uint64_t svlen = Nmax / (30ull * 8) + 1;

    ParallelSieve(Nmax, vPrimes, svlen);

    return CountAllSegments(vPrimes, Nmax, svlen);
    //return CountAllOnce(vPrimes, Nmax, svlen);
    //return 11;
}
