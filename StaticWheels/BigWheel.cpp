#include "Helper.h"
#include "wheels.h"


// *****************************************************
//  Big Wheel
// *****************************************************

wheels::BigWheel7 bigwheel;
constexpr uint64_t chunk_span = chunk_bits * bigwheel.wheel_span;

void testbigdiv(void)
{
    uint64_t t = 0;
    for (unsigned i = 0; i < 0xffff'ffff; i++)
        t += i / bigwheel.wheel_span;
    std::cout << t;
}

extern uint64_t segments_no, segment_bytes, segment_span;

unsigned first_root_prime_idx7;
unsigned first_root_prime7 = bigwheel.wheel[1];
extern uint8_t root_primes_gap[];

unsigned rootsieveR7[NO_ROOT_PRIMES];

void ResetSeed(void);
unsigned NextBatchT(void);
unsigned NextBatchP(void);

extern thread_local uint64_t* vlastpk7 = NULL;
// local sieve
extern thread_local uint8_t* sv7 = NULL;

void SieveChunk7(const uint64_t seg_start, const uint64_t stop, const uint64_t val)
{
    const uint64_t maxn = std::min(segment_bytes * byte_bits, (stop - seg_start) / bigwheel.wheel_span + 1);
    const uint64_t maxp = (uint64_t)floor(sqrt(stop));
    uint64_t p = first_root_prime7;
    unsigned i = first_root_prime_idx7;
    //for (unsigned i = first_root_prime_idx7; i < NO_ROOT_PRIMES;)
    while(p <= maxp)
    {
        tpPrime n = vlastpk7[i];
        //if (n < maxn)
        {
            //uint64_t idx0 = n / 8, bit0 = n % 8;
            //uint64_t idxdeltap = p / 8, idxdeltab = p % 8;
            for (; n < maxn; n += p)
            {
                //uint64_t idx = n / 8, bit = n % 8;
                //uint64_t nn = seg_start + wheel.wheel_span * (8 * idx + bit) + val;
                sv7[n / byte_bits] |= BIT_MASK[n % byte_bits];

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
            vlastpk7[i] = n;
        }
        p += root_primes_gap[++i]; //if (p > maxp) break;
    }
}

void InitLastPK7(const uint64_t start, const uint64_t st, const unsigned ridx)
{
    const uint64_t maxp = (uint64_t)floor(sqrt(st));
    uint64_t p = first_root_prime7;
    unsigned i = first_root_prime_idx7;
    while (p <= maxp) // (i < NO_ROOT_PRIMES)
    {   // for each root prime
        // get initial multiplier approximation
        const uint64_t startk = (start > p * p) ? (start / p) : p;
        // adjust to multiple of wheel_span
        uint64_t k = bigwheel.wheel_span * (startk / bigwheel.wheel_span);
        // move to current factor
        k += (*bigwheel.factor[ridx])[rootsieveR7[i]];
        k += bigwheel.wheel_span * (k < startk); //if (k < startk) k += bigwheel.wheel_span;
        uint64_t n = k * p;
        // adjust if lower
        while (n < start) { k += bigwheel.wheel_span; n = k * p; }
        vlastpk7[i] = (n - start) / bigwheel.wheel_span;
        // get next root prime value
        p += root_primes_gap[++i]; //if (p > maxp) break;
    }
    assert(i < NO_ROOT_PRIMES); // should practically never reach the end of the list
}

unsigned root_k[NO_ROOT_PRIMES + 1];
void InitLastPK7x(/*const uint64_t start, */const uint64_t st, const unsigned ridx)
{
    //assert(start == 0);
    const uint64_t maxp = (uint64_t)floor(sqrt(st));
    uint64_t p = first_root_prime7;
    unsigned i = first_root_prime_idx7;
    while (p <= maxp) // (i < NO_ROOT_PRIMES)
    {   // for each root prime
        // get initial multiplier approximation
        uint64_t k = root_k[i] + (*bigwheel.factor[ridx])[rootsieveR7[i]];
        k += bigwheel.wheel_span * (k < p); 
        vlastpk7[i] = (k * p) / bigwheel.wheel_span;
        // get next root prime value
        p += root_primes_gap[++i]; //if (p > maxp) break;
    }
    assert(i < NO_ROOT_PRIMES); // should practically never reach the end of the list
}

uint64_t ProcessBatch7(const uint64_t Nmax, const unsigned batch)
{
    const uint64_t segment = batch / bigwheel.wheel_values;
    const unsigned ridx = batch % bigwheel.wheel_values;
    const uint64_t seg_start = segment * segment_span; if (seg_start > Nmax) return 0;
    const uint64_t seg_stop = std::min(Nmax, seg_start + segment_span);
    const unsigned value = bigwheel.wheel[ridx];

    //InitLastPK7(seg_start, seg_stop, ridx);
    InitLastPK7x(/*seg_start, */seg_stop, ridx);
    memset(sv7, 0, segment_bytes);
    //std::cout << imax << " | " << segment_bytes << " | " << segment_span; nln();

    for (uint64_t ck_start = seg_start; ck_start < seg_stop; ck_start += chunk_span)
    {
        const uint64_t ck_stop = std::min(ck_start + chunk_span, seg_stop);
        SieveChunk7(seg_start, ck_stop, value);
    }

    uint64_t numPrimes = 0;
    const unsigned imax = unsigned((seg_stop - seg_start) / byte_bits / bigwheel.wheel_span);
    for (unsigned i = 0; i < imax; i++)
        numPrimes += bit0_count_table[sv7[i]];
    if (imax < segment_bytes)
        for (unsigned b = 0; b < byte_bits; b++)
            if (!(sv7[imax] & BIT_MASK[b]))
            {
                uint64_t n = seg_start + bigwheel.wheel_span * ((uint64_t)byte_bits * imax + b) + value;
                if (n > Nmax) break;
                numPrimes++;
            }
    //std::cout << seg_start << " | " << seg_stop << " | " << value; nln();
    return numPrimes;
}

uint64_t ProcessThread7(const uint64_t Nmax)
{
    // for each segment|value pair we have one batch 
    assert(segments_no * bigwheel.wheel_values < UINT32_MAX);
    const unsigned batches_no = segments_no * bigwheel.wheel_values;

    vlastpk7 = new tpPrime[NO_ROOT_PRIMES]; //vlastpk[0] = 0;
    sv7 = new uint8_t[segment_bytes];
    //std::cout << std::hex << (uint64_t)(sieve) << " | "; nln();

    uint64_t numPrimes = 0;
    for (unsigned batch = NextBatchT(); batch < batches_no; batch = NextBatchT())
    { // do one batch
        numPrimes += ProcessBatch7(Nmax, batch);
#define SHOW_PROGRESS
#ifdef SHOW_PROGRESS
        if(batch % 32 == 0)
            fprintf(stdout, "%02u%%\b\b\b", batch * 100 / batches_no);
            //fprintf(stdout, "%5.1f%%\b\b\b\b\b\b", batch * 100.0 / batches_no);
#endif
    }
#ifdef SHOW_PROGRESS
        fprintf(stdout, "        ");
#endif
    delete[] vlastpk7;
    delete[] sv7;

    return numPrimes;
}

void SetSegmentsNo7(const uint64_t Nmax, const unsigned nativeThreads)
{
    // total chunks
    uint64_t chunks = Nmax / chunk_span;
    if (Nmax > chunks * chunk_span) chunks++;
    // we have at least one thread per value 
    // and we strive to spread batches on cores 
    unsigned l = std::lcm(nativeThreads, bigwheel.wheel_values);
    unsigned segments = l / bigwheel.wheel_values;
    uint64_t segment_chuncks = chunks / segments;
    if (chunks > segments * segment_chuncks) segment_chuncks++;
    segment_bytes = segment_chuncks * chunk_bytes;
    while (segment_bytes > max_segment_bytes)  // we need to adjust
    {
        // get next value that spreads equally
        for (segments++; segments * bigwheel.wheel_values % nativeThreads > 0; segments++);
        segment_chuncks = chunks / segments;
        if (chunks > segments * segment_chuncks) segment_chuncks++;
        segment_bytes = segment_chuncks * chunk_bytes;
    }
    segments_no = segments;
    segment_span = segment_bytes * bigwheel.wheel_span * byte_bits;
    std::cout << "[ " << nativeThreads << " threads; " << segments_no;
    std::cout << " segment(s); " << bigwheel.wheel_values << " values ]";
}

uint64_t ParallelSieve7(const uint64_t Nmax)
{
    const unsigned nativeThreads = 1 * std::max(1u, std::thread::hardware_concurrency());

    ResetSeed();
    SetSegmentsNo7(Nmax, nativeThreads);

    uint64_t numPrimes = first_root_prime_idx7 + 1ull;
    std::vector<std::future<uint64_t>> counters;
    for (unsigned t = 0; t < nativeThreads; t++)
    {
        counters.push_back(std::async(&ProcessThread7, Nmax));
        //numPrimes += ProcessThread7(Nmax);
    }

    for (auto& c : counters)
        numPrimes += c.get();

    return numPrimes;
}

void SoA_LP_gen_root_primes(void);

void InitRootPrimes7(void)
{
    SoA_LP_gen_root_primes();
    unsigned p = 5, i = 1;
    // skip first root primes excluded by the wheel
    while (p < first_root_prime7)
#if LIMIT_lg10 < 9
        p += root_primes_gap[i++];
#else
        p += root_primes_gap[i++] * 2;
#endif
    first_root_prime_idx7 = --i;
    while (i < NO_ROOT_PRIMES)
    {
        const unsigned r = p % bigwheel.wheel_span;
        if (r == 1)
            rootsieveR7[i] = 0;
        else
        {
            auto itr = std::lower_bound(bigwheel.wheel.data() + 1, bigwheel.wheel.data() + bigwheel.wheel_values, r);
            rootsieveR7[i] = (unsigned)(itr - bigwheel.wheel.data());
            assert((itr - bigwheel.wheel.data()) < std::numeric_limits<unsigned>::max());
            assert(bigwheel.wheel[rootsieveR7[i]] == r);
        }
        root_k[i] = bigwheel.wheel_span * (p / bigwheel.wheel_span);
#if LIMIT_lg10 < 9
        p += root_primes_gap[++i];
#else
        p += root_primes_gap[++i] * 2;
#endif
    }
}
uint64_t SoE_pat_w7(const uint64_t Nmax, void*, void*, void*)
{
    InitRootPrimes7();
    return ParallelSieve7(Nmax);
}

void InitBigWheel(void)
{
    nln(true); bigwheel.GenFactors();
}



constexpr auto SZ = bigwheel.wheel_values;
constexpr auto BLKSZ = 256u; static_assert(SZ % BLKSZ == 0);
std::array<unsigned*, SZ> arr;
void testalloc(void)
{
    cTimer tmr; tmr.Start();

    std::for_each(std::execution::par_unseq,
        range(SZ / BLKSZ).begin(), range(SZ / BLKSZ).end(), [&](const unsigned i)
        {
            unsigned *pv = new unsigned[BLKSZ * SZ];
            for (unsigned j = 0; j < BLKSZ; j++) 
                arr[i*BLKSZ + j] = pv + j;
        });

    std::for_each(std::execution::par_unseq,
        range(SZ / BLKSZ).begin(), range(SZ / BLKSZ).end(), [&](const unsigned i)
        {
            delete arr[i * BLKSZ];
        });

    tmr.Stop(true, "mem");
}

