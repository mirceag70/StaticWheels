#include "Helper.h"
#include "wheels.h"

constexpr unsigned WHEEL_N = 6; // [3..6]
inline wheels::RootSieve<WHEEL_N> wheel;
constexpr uint64_t chunk_span = chunk_bits * wheel.wheel_span;

void test_wheels(void)
{
    cTimer tmr; tmr.Start();
    //for (auto j : range(100))
    //{
    //    wheel.GenMods();
    //}
    nln();
    for (auto i : range(wheel.wheel_values))
    {
        std::cout << i << "\n";
        for (auto j : range(wheel.wheel_values))
            std::cout << +wheel.mods[i][j] << " | ";
        nln();
        for (auto j : range(wheel.wheel_values))
            std::cout << +wheel.factor[i][j] << " | ";
        nln();
    }
    tmr.Stop(true, "mods");
}

void SoA_LP_gen_root_primes(void);
extern uint8_t root_primes_gap[];
extern uint8_t rsvR[];
typedef wheels::RootSieve<WHEEL_N>::tpWheelEntry rootsieveRtp;
rootsieveRtp rootsieveR[NO_ROOT_PRIMES];
unsigned first_root_prime = wheel.wheel[1];
unsigned first_root_prime_idx;
unsigned szi1, szi3, szi5, szi7;
void InitRootPrimes(void)
{
    SoA_LP_gen_root_primes();
    unsigned p = 5, i = 1;
    szi1 = szi3 = szi5 = szi7 = 0;
    // skip first root primes excluded by the wheel
    while (p < first_root_prime)
#if LIMIT_lg10 < 9
        p += root_primes_gap[i++];
#else
        p += root_primes_gap[i++] * 2;
#endif
    first_root_prime_idx = --i;
    while(i < NO_ROOT_PRIMES)
    {
        //switch (p % 8) { default: NEVERHERE;
        //case 1: szi1++; break; case 3: szi3++; break;
        //case 5: szi5++; break; case 7: szi7++; break; };

        const unsigned r = p % wheel.wheel_span;
        if (r == 1)
            rootsieveR[i] = 0;
        else
        {
            auto itr = std::lower_bound(wheel.wheel.data() + 1, wheel.wheel.data() + wheel.wheel_values, r);
            rootsieveR[i] = (rootsieveRtp)(itr - wheel.wheel.data());
            assert((itr - wheel.wheel.data()) < std::numeric_limits<rootsieveRtp>::max());
            assert(wheel.wheel[rootsieveR[i]] == r);
            //assert(rsvR[i] == rootsieveR[i]);
        }
#if LIMIT_lg10 < 9
        p += root_primes_gap[++i];
#else
        p += root_primes_gap[++i] * 2;
#endif
    }
}

// vector with last multiples of root primes p * k
thread_local uint64_t* vlastpk = NULL;
// local sieve
thread_local uint8_t* sieve = NULL;

void InitLastPK(const uint64_t start, const uint64_t stop, const unsigned ridx)
{
    const uint64_t maxp = (uint64_t)floor(sqrt(stop));
    uint64_t p = first_root_prime;
    unsigned i = first_root_prime_idx;
    while (i < NO_ROOT_PRIMES)
    {   // for each root prime
        // get initial multiplier approximation
        const uint64_t startk = (start > p * p) ? (start / p) : p;
        // adjust to multiple of wheel_span
        uint64_t k = wheel.wheel_span * (startk / wheel.wheel_span);
        // move to current factor
        k += wheel.factor[ridx][rootsieveR[i]];
        //k += wheel.factor[rootsieveR[i]][ridx];
        k += wheel.wheel_span * (k < startk); //if (k < startk) k += wheel.wheel_span;
        uint64_t n = k * p;
        // adjust if lower
        while (n < start) { k += wheel.wheel_span; n = k * p; }
        vlastpk[i] = (n - start) / wheel.wheel_span;
        // get next root prime value
        p += root_primes_gap[++i]; if (p > maxp) break;
    }
    assert(i < NO_ROOT_PRIMES); // should practically never reach the end of the list
}

//__forceinline
//__declspec(noinline)
void InitLastPKxx(const uint64_t stop, const unsigned ridx)
{
    //if (start > 0) __debugbreak();
    const uint64_t maxp = (uint64_t)floor(sqrt(stop));
    uint64_t p = first_root_prime;
    unsigned i = first_root_prime_idx;
    while (i < NO_ROOT_PRIMES)
    {   // for each root prime
        // get initial multiplier approximation
        //const uint64_t startk = /*(start > p * p) ? (start / p) :*/ p;
        //if (startk != p) __debugbreak();
        // adjust to multiple of wheel_span
        //uint64_t k = 0;// wheel.wheel_span* (startk / wheel.wheel_span);
        // move to current factor
        uint64_t k = wheel.factor[ridx][rootsieveR[i]];
        //k += wheel.factor[rootsieveR[i]][ridx];
        k += wheel.wheel_span * (k < p); //if (k < startk) k += wheel.wheel_span;
        //uint64_t n = k * p;
        // adjust if lower
        //while (n < start) { k += wheel.wheel_span; n = k * p; }
        vlastpk[i] = k * p / wheel.wheel_span;
        // get next root prime value
        p += root_primes_gap[++i]; if (p > maxp) break;
    }
    assert(i < NO_ROOT_PRIMES); // should practically never reach the end of the list
}

uint64_t segments_no, segment_bytes, segment_span;

void SieveChunk(const uint64_t seg_start, const uint64_t stop, const uint64_t val)
{
    const uint64_t maxn = std::min (segment_bytes * byte_bits, (stop - seg_start) / wheel.wheel_span + 1);
    const uint64_t maxp = (uint64_t)floor(sqrt(stop));
    uint64_t p = first_root_prime;
    unsigned i = first_root_prime_idx;
    //for (unsigned i = first_root_prime_idx; i < NO_ROOT_PRIMES;)
    while(p <= maxp)
    {
        tpPrime n = vlastpk[i];
        //if (n < maxn)
        {
            //uint64_t idx0 = n / 8, bit0 = n % 8;
            //uint64_t idxdeltap = p / 8, idxdeltab = p % 8;
            for (; n < maxn; n += p)
            {
                //uint64_t idx = n / 8, bit = n % 8;
                //uint64_t nn = seg_start + wheel.wheel_span * (8 * idx + bit) + val;
                sieve[n / byte_bits] |= BIT_MASK[n % byte_bits];
                
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
            vlastpk[i] = n;
        }
        p += root_primes_gap[++i]; //if (p > maxp) break;
    }
}

uint64_t ProcessBatch(const uint64_t Nmax, const unsigned batch)
{
    const uint64_t segment = batch / wheel.wheel_values;
    const unsigned ridx = batch % wheel.wheel_values;
    const uint64_t seg_start = segment * segment_span; if (seg_start > Nmax) return 0;
    const uint64_t seg_stop = std::min(Nmax, seg_start + segment_span);
    const unsigned value = wheel.wheel[ridx];

    InitLastPK(seg_start, seg_stop, ridx);
    memset(sieve, 0, segment_bytes);
    //std::cout << imax << " | " << segment_bytes << " | " << segment_span; nln();

    for (uint64_t ck_start = seg_start; ck_start < seg_stop; ck_start += chunk_span)
    {
        const uint64_t ck_stop = std::min(ck_start + chunk_span, seg_stop);
        SieveChunk(seg_start, ck_stop, value);
    }

    uint64_t numPrimes = 0;
    const unsigned imax = unsigned((seg_stop - seg_start) / byte_bits / wheel.wheel_span);
    for (unsigned i = 0; i < imax; i++)
        numPrimes += bit0_count_table[sieve[i]];
    if(imax < segment_bytes)
    for (unsigned b = 0; b < byte_bits; b++)
        if (!(sieve[imax] & BIT_MASK[b]))
        {
            uint64_t n = seg_start + wheel.wheel_span * ((uint64_t)byte_bits * imax + b) + value;
            if (n > Nmax) break;
            numPrimes++;
        }
    //std::cout << seg_start << " | " << seg_stop << " | " << value; nln();
    return numPrimes;
}


std::atomic<unsigned> tseed;    // batch seed - for thread id
std::atomic<unsigned> pseed;    // batch seed - for progress
void ResetSeed(void) { tseed.store(0); pseed.store(0); }
unsigned NextBatchT(void) { return tseed.fetch_add(1u, std::memory_order_relaxed); }
unsigned NextBatchP(void) { return pseed.fetch_add(1u, std::memory_order_relaxed); }

uint64_t ProcessThread(const uint64_t Nmax)
{
    // for each segment|value pair we have one batch 
    assert(segments_no * bigwheel.wheel_values < UINT32_MAX);
    const unsigned batches_no = segments_no * wheel.wheel_values;

    vlastpk = new tpPrime[NO_ROOT_PRIMES]; vlastpk[0] = 0;
    sieve = new uint8_t[segment_bytes];
    //std::cout << std::hex << (uint64_t)(sieve) << " | "; nln();

    uint64_t numPrimes = 0;
    for (unsigned batch = NextBatchT(); batch < batches_no; batch = NextBatchT())
    { // do one batch
        numPrimes += ProcessBatch(Nmax, batch);
#define SHOW_PROGRESS
#ifdef SHOW_PROGRESS
        if (batch % 32 == 0)
            fprintf(stdout, "%2u%%\b\b\b", batch * 100 / batches_no);
            //fprintf(stdout, "%5.1f%%\b\b\b\b\b\b", batch * 100.0 / batches_no);
#endif        
    }
#ifdef SHOW_PROGRESS
    fprintf(stdout, "        ");
#endif'
    delete[] vlastpk;
    delete[] sieve;

    return numPrimes;
}

void SetSegmentsNo(const uint64_t Nmax, const unsigned nativeThreads)
{
    // total chunks
    uint64_t chunks = Nmax / chunk_span;
    if (Nmax > chunks * chunk_span) chunks++;
    // we have at least one thread per value 
    // and we strive to spread batches on cores 
    unsigned l = std::lcm(nativeThreads, wheel.wheel_values);
    unsigned segments = l / wheel.wheel_values;
    uint64_t segment_chuncks = chunks / segments;
    if (chunks > segments * segment_chuncks) segment_chuncks++;
    segment_bytes = segment_chuncks * chunk_bytes;
    while (segment_bytes > max_segment_bytes)  // we need to adjust
    {   
        // get next value that spreads equally
        for (segments++; segments * wheel.wheel_values % nativeThreads > 0; segments++);
        segment_chuncks = chunks / segments;
        if (chunks > segments * segment_chuncks) segment_chuncks++;
        segment_bytes = segment_chuncks * chunk_bytes;
    }
    segments_no = segments;
    segment_span = segment_bytes * wheel.wheel_span * byte_bits;
    std::cout << "[ " << nativeThreads << " threads; " << segments_no;
    std::cout << " segment(s); W" << WHEEL_N << "-" << wheel.wheel_values << " values ]";
}

uint64_t ParallelSieve(const uint64_t Nmax)
{
    const unsigned nativeThreads = 1*std::max(1u, std::thread::hardware_concurrency());

    ResetSeed();
    SetSegmentsNo(Nmax, nativeThreads);

    uint64_t numPrimes = first_root_prime_idx + 1ull;
    std::vector<std::future<uint64_t>> counters;
    for (unsigned t = 0; t < nativeThreads; t++)
    {
        counters.push_back(std::async(&ProcessThread, Nmax));
        //numPrimes += ProcessThread(Nmax);
    }

    for (auto& c : counters) 
        numPrimes += c.get();

    return numPrimes;
}

uint64_t SoE_pat_free(const uint64_t Nmax, void*, void*, void*)
{
    InitRootPrimes();
    return ParallelSieve(Nmax);
}



// *****************************************************
//  Bit reset pattern
// *****************************************************

void SieveChunkBx0(const uint64_t seg_start, const uint64_t stop, const uint64_t val)
{
    const uint64_t maxn = std::min(segment_bytes * byte_bits, (stop - seg_start) / wheel.wheel_span + 1);
    const uint64_t maxp = (uint64_t)floor(sqrt(stop));
    uint64_t p = first_root_prime;
    for (unsigned i = first_root_prime_idx; i < NO_ROOT_PRIMES;)
    {
        tpPrime n = vlastpk[i];
        if (n < maxn)
        {
            uint64_t idx0 = n / 8, bit0 = n % 8;
            uint64_t idxdeltap = p / 8, idxdeltab = p % 8;
            for (; n < maxn; n += p)
            {
                //uint64_t idx = n / 8, bit = n % 8;
                //uint64_t nn = seg_start + wheel.wheel_span * (8 * idx + bit) + val;
                sieve[n / byte_bits] |= BIT_MASK[n % byte_bits];

                //sv[idx] |= BIT_MASK[bit];
                uint64_t idx = n / 8, bit = n % 8;
                assert(idx0 == idx and bit0 == bit);
                //sieve[idx0] |= BIT_MASK[bit0];
                switch (idxdeltab)
                {
                default: NEVERHERE;
                case 1: 
                    if (bit0 == 7) { bit0 = 0; idx0++; } 
                    else bit0++; break;
                case 3: 
                    if (bit0 < 5) { bit0 += 3; }
                    else { bit0-=5; idx0++; }; break;
                case 5:
                    if (bit0 < 3) { bit0 += 5; }
                    else { bit0 -= 3; idx0++; }; break;
                case 7:
                    if (bit0 == 0) { bit0 = 7; }
                    else { bit0--; idx0++; }; break;
                }
                idx0 += idxdeltap;
            }
            vlastpk[i] = n;
        }
        p += root_primes_gap[++i]; if (p > maxp) break;
    }
}

void SieveChunkB(const uint64_t seg_start, const uint64_t stop, const uint64_t val)
{
    const uint64_t maxn = std::min(segment_bytes * byte_bits, (stop - seg_start) / wheel.wheel_span + 1);
    const uint64_t maxp = (uint64_t)floor(sqrt(stop));
    uint64_t p = first_root_prime;
    for (unsigned i = first_root_prime_idx; i < NO_ROOT_PRIMES;)
    {
        tpPrime n = vlastpk[i];
        if (n < maxn)
        {
            uint64_t idx0 = n / 8, bit0 = n % 8;
            uint64_t idxdeltap = p / 8, idxdeltab = p % 8;
            /*for (; n < maxn; n += p)
            {
                //sieve[n / byte_bits] |= BIT_MASK[n % byte_bits];

                //uint64_t idx = n / 8, bit = n % 8;
                //assert(idx0 == idx and bit0 == bit);
                //uint64_t nn = seg_start + wheel.wheel_span * (8 * idx + bit) + val;
                sieve[idx0] |= BIT_MASK[bit0];
                switch (idxdeltab)
                {
                default: NEVERHERE;
                case 1:
                    if (bit0 == 7) { bit0 = 0; idx0++; }
                    else bit0++; break;
                case 3:
                    if (bit0 < 5) { bit0 += 3; }
                    else { bit0 -= 5; idx0++; }; break;
                case 5:
                    if (bit0 < 3) { bit0 += 5; }
                    else { bit0 -= 3; idx0++; }; break;
                case 7:
                    if (bit0 == 0) { bit0 = 7; }
                    else { bit0--; idx0++; }; break;
                }
                idx0 += idxdeltap;
            }*/
            switch (idxdeltab)
            {
            default: NEVERHERE;
            case 1:
                //for (; n < maxn; n += p)
                //{
                //    if (n >= maxn) goto SWITCH_OUT;
                //    sieve[idx0] |= BIT_MASK[bit0];
                //    if (bit0 == 7) { bit0 = 0; idx0++; }
                //    else bit0++;
                //    idx0 += idxdeltap;
                //}; break;
                switch (bit0)
                {
                default: NEVERHERE;
                    for (;;)
                    {
                case 0:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[0]; idx0 += idxdeltap;
                case 1:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[1]; idx0 += idxdeltap;
                case 2:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[2]; idx0 += idxdeltap;
                case 3:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[3]; idx0 += idxdeltap;
                case 4:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[4]; idx0 += idxdeltap;
                case 5:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[5]; idx0 += idxdeltap;
                case 6:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[6]; idx0 += idxdeltap;
                case 7:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[7]; idx0 += idxdeltap + 1;
                    };
                } break;
            case 3:
                //for (; n < maxn; n += p)
                //{
                //    sieve[idx0] |= BIT_MASK[bit0];
                //    if (bit0 < 5) { bit0 += 3; }
                //    else { bit0 -= 5; idx0++; };
                //    idx0 += idxdeltap;
                //}; break;
                switch (bit0)
                {
                default: NEVERHERE;
                    for (;;)
                    {
                case 0:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[0]; idx0 += idxdeltap;
                case 3:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[3]; idx0 += idxdeltap;
                case 6:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[6]; idx0 += idxdeltap + 1;
                case 1:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[1]; idx0 += idxdeltap;
                case 4:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[4]; idx0 += idxdeltap;
                case 7:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[7]; idx0 += idxdeltap + 1;
                case 2:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[2]; idx0 += idxdeltap;
                case 5:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[5]; idx0 += idxdeltap + 1;
                    }
                } break;
            case 5:
                //for (; n < maxn; n += p)
                //{
                //    sieve[idx0] |= BIT_MASK[bit0];
                //    if (bit0 < 3) { bit0 += 5; }
                //    else { bit0 -= 3; idx0++; };
                //    idx0 += idxdeltap;
                //};  break;
                switch (bit0)
                {
                default: NEVERHERE;
                    for (;;)
                    {
                case 0:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[0]; idx0 += idxdeltap;
                case 5:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[5]; idx0 += idxdeltap + 1;
                case 2:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[2]; idx0 += idxdeltap;
                case 7:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[7]; idx0 += idxdeltap + 1;
                case 4:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[4]; idx0 += idxdeltap + 1;
                case 1:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[1]; idx0 += idxdeltap;
                case 6:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[6]; idx0 += idxdeltap + 1;
                case 3:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[3]; idx0 += idxdeltap + 1;
                    }
                } break;
            case 7:
                //for (; n < maxn; n += p)
                //{
                //    sieve[idx0] |= BIT_MASK[bit0];
                //    if (bit0 == 0) { bit0 = 7; }
                //    else { bit0--; idx0++; };
                //    idx0 += idxdeltap;
                //};  break;
                switch (bit0)
                {
                default: NEVERHERE;
                    for (;;)
                    {
                case 0:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[0]; idx0 += idxdeltap;
                case 7:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[7]; idx0 += idxdeltap + 1;
                case 6:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[6]; idx0 += idxdeltap + 1;
                case 5:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[5]; idx0 += idxdeltap + 1;
                case 4:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[4]; idx0 += idxdeltap + 1;
                case 3:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[3]; idx0 += idxdeltap + 1;
                case 2:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[2]; idx0 += idxdeltap + 1;
                case 1:
                    if (n >= maxn) goto SWITCH_OUT; else n += p;
                    sieve[idx0] |= BIT_MASK[1]; idx0 += idxdeltap + 1;
                    };
                } break;
            };
            SWITCH_OUT: vlastpk[i] = n;
        }
        p += root_primes_gap[++i]; if (p > maxp) break;
    }
}

void InitLastPKB(const uint64_t start, const uint64_t stop, const unsigned ridx)
{
    unsigned sml, med, big;
    sml = med = big = 0;
    const uint64_t maxp = (uint64_t)floor(sqrt(stop));
    uint64_t p = first_root_prime;
    unsigned i = first_root_prime_idx;
    while (i < NO_ROOT_PRIMES)
    {   // for each root prime
        // get initial multiplier approximation
        const uint64_t startk = (start > p * p) ? (start / p) : p;
        // adjust to multiple of wheel_span
        uint64_t k = wheel.wheel_span * (startk / wheel.wheel_span);
        if (p < chunk_span)     sml++; else
        if (p < segment_span)   med++; else
                                big++;
        // move to current factor
        k += wheel.factor[rootsieveR[i]][ridx];
        if (k < startk) k += wheel.wheel_span;
        uint64_t n = k * p;
        // adjust if lower
        while (n < start) { k += wheel.wheel_span; n = k * p; }
        vlastpk[i] = (n - start) / wheel.wheel_span;
        // get next root prime value
        p += root_primes_gap[++i]; if (p > maxp) break;
    }
    assert(i < NO_ROOT_PRIMES); // should practically never reach the end of the list
    assert(med == 0 and big == 0);
}

uint64_t ProcessBatchB(const uint64_t Nmax, const unsigned batch)
{
    const uint64_t segment = batch / wheel.wheel_values;
    const unsigned ridx = batch % wheel.wheel_values;
    const uint64_t seg_start = segment * segment_span; if (seg_start > Nmax) return 0;
    const uint64_t seg_stop = std::min(Nmax, seg_start + segment_span);
    const unsigned value = wheel.wheel[ridx];

    InitLastPKB(seg_start, seg_stop, ridx);
    memset(sieve, 0, segment_bytes);
    //std::cout << imax << " | " << segment_bytes << " | " << segment_span; nln();

    for (uint64_t ck_start = seg_start; ck_start < seg_stop; ck_start += chunk_span)
    {
        const uint64_t ck_stop = std::min(ck_start + chunk_span, seg_stop);
        SieveChunkB(seg_start, ck_stop, value);
    }

    uint64_t numPrimes = 0;
    const unsigned imax = unsigned((seg_stop - seg_start) / byte_bits / wheel.wheel_span);
    for (unsigned i = 0; i < imax; i++)
        numPrimes += bit0_count_table[sieve[i]];
    if (imax < segment_bytes)
        for (unsigned b = 0; b < byte_bits; b++)
            if (!(sieve[imax] & BIT_MASK[b]))
            {
                uint64_t n = seg_start + wheel.wheel_span * ((uint64_t)byte_bits * imax + b) + value;
                if (n > Nmax) break;
                numPrimes++;
            }
    //std::cout << seg_start << " | " << seg_stop << " | " << value; nln();
    return numPrimes;
}

uint64_t ProcessThreadB(const uint64_t Nmax)
{
    // for each segment|value pair we have one batch 
    const uint64_t batches_no = segments_no * wheel.wheel_values;

    vlastpk = new tpPrime[NO_ROOT_PRIMES]; vlastpk[0] = 0;
    sieve = new uint8_t[segment_bytes];
    //std::cout << std::hex << (uint64_t)(sieve) << " | "; nln();

    uint64_t numPrimes = 0;
    for (unsigned batch = NextBatchT(); batch < batches_no; batch = NextBatchT())
    { // do one batch
        numPrimes += ProcessBatchB(Nmax, batch);
        fprintf(stdout, "%5.1f%%\b\b\b\b\b\b", NextBatchP() * 100.0 / batches_no);
    }
    delete[] vlastpk;
    delete[] sieve;

    return numPrimes;
}

uint64_t ParallelSieveB(const uint64_t Nmax)
{
    const unsigned nativeThreads = std::max(1u, std::thread::hardware_concurrency());

    ResetSeed();
    SetSegmentsNo(Nmax, nativeThreads);

    uint64_t numPrimes = first_root_prime_idx + 1ull;
    std::vector<std::future<uint64_t>> counters;
    for (unsigned t = 0; t < nativeThreads; t++)
    {
        counters.push_back(std::async(&ProcessThreadB, Nmax));
        //numPrimes += ProcessThreadB(Nmax);
    }

    for (auto& c : counters)
        numPrimes += c.get();

    return numPrimes;
}

uint64_t SoE_bitpat_free0(const uint64_t Nmax, void*, void*, void*)
{
    InitRootPrimes();
    return ParallelSieveB(Nmax);
}

// *****************************************************
//  Bit reset pattern + list init
// *****************************************************

struct /*alignas(64)*/ lastpk
{
    unsigned idx;
    unsigned delta;
    short unsigned bit;
    short unsigned _pad;
};

thread_local unsigned maxi1, maxi3, maxi5, maxi7;
thread_local lastpk * vlastpk1, * vlastpk3, * vlastpk5, * vlastpk7;

void SieveChunkB1(const uint64_t seg_start, const uint64_t stop, const uint64_t val)
{
    const uint64_t maxi = std::min(segment_bytes, (stop - seg_start) / byte_bits / wheel.wheel_span + 1);
    const uint64_t maxp = (uint64_t)floor(sqrt(stop / byte_bits));
    for (unsigned i = 0; i < maxi1; i++)
    {
        unsigned idx0 = vlastpk1[i].idx, delta0 = vlastpk1[i].delta, delta1 = delta0 + 1;
        //if (vlastpk1[i].bit > 3)
        //    if (vlastpk1[i].bit > 5)
        //        if (vlastpk1[i].bit == 6) goto CASE_1_6;
        //        else { assert(vlastpk1[i].bit == 7); goto CASE_1_7; }
        //    else
        //        if (vlastpk1[i].bit == 4) goto CASE_1_4;
        //        else { assert(vlastpk1[i].bit == 5); goto CASE_1_5; }
        //else
        //    if (vlastpk1[i].bit > 1)
        //        if (vlastpk1[i].bit == 2) goto CASE_1_2;
        //        else { assert(vlastpk1[i].bit == 3); goto CASE_1_3; }
        //    else
        //        if (vlastpk1[i].bit == 0) goto CASE_1_0;
        //        else { assert(vlastpk1[i].bit == 1); goto CASE_1_1; }
        ////loop
        //CASE_1_0:
        //    if (idx0 >= maxi) { vlastpk1[i].bit = 0; goto SWITCH_OUT1; }
        //    sieve[idx0] |= BIT_MASK[0]; idx0 += delta0;
        //CASE_1_1:
        //    if (idx0 >= maxi) { vlastpk1[i].bit = 1; goto SWITCH_OUT1; }
        //    sieve[idx0] |= BIT_MASK[1]; idx0 += delta0;
        //CASE_1_2:
        //    if (idx0 >= maxi) { vlastpk1[i].bit = 2; goto SWITCH_OUT1; }
        //    sieve[idx0] |= BIT_MASK[2]; idx0 += delta0;
        //CASE_1_3:
        //    if (idx0 >= maxi) { vlastpk1[i].bit = 3; goto SWITCH_OUT1; }
        //    sieve[idx0] |= BIT_MASK[3]; idx0 += delta0;
        //CASE_1_4:
        //    if (idx0 >= maxi) { vlastpk1[i].bit = 4; goto SWITCH_OUT1; }
        //    sieve[idx0] |= BIT_MASK[4]; idx0 += delta0;
        //CASE_1_5:
        //    if (idx0 >= maxi) { vlastpk1[i].bit = 5; goto SWITCH_OUT1; }
        //    sieve[idx0] |= BIT_MASK[5]; idx0 += delta0;
        //CASE_1_6:
        //    if (idx0 >= maxi) { vlastpk1[i].bit = 6; goto SWITCH_OUT1; }
        //    sieve[idx0] |= BIT_MASK[6]; idx0 += delta0;
        //CASE_1_7:
        //    if (idx0 >= maxi) { vlastpk1[i].bit = 7; goto SWITCH_OUT1; }
        //    sieve[idx0] |= BIT_MASK[7]; idx0 += delta0 + 1;
        ////loop back
        //goto CASE_1_0;
        switch (vlastpk1[i].bit)
        {
            default: NEVERHERE;
                for (;;)
                {
            case 0:
                if (idx0 >= maxi) { vlastpk1[i].bit = 0; goto SWITCH_OUT1; }
                sieve[idx0] |= BIT_MASK[0]; idx0 += delta0;
            case 1:
                if (idx0 >= maxi) { vlastpk1[i].bit = 1; goto SWITCH_OUT1; }
                sieve[idx0] |= BIT_MASK[1]; idx0 += delta0;
            case 2:
                if (idx0 >= maxi) { vlastpk1[i].bit = 2; goto SWITCH_OUT1; }
                sieve[idx0] |= BIT_MASK[2]; idx0 += delta0;
            case 3:
                if (idx0 >= maxi) { vlastpk1[i].bit = 3; goto SWITCH_OUT1; }
                sieve[idx0] |= BIT_MASK[3]; idx0 += delta0;
            case 4:
                if (idx0 >= maxi) { vlastpk1[i].bit = 4; goto SWITCH_OUT1; }
                sieve[idx0] |= BIT_MASK[4]; idx0 += delta0;
            case 5:
                if (idx0 >= maxi) { vlastpk1[i].bit = 5; goto SWITCH_OUT1; }
                sieve[idx0] |= BIT_MASK[5]; idx0 += delta0;
            case 6:
                if (idx0 >= maxi) { vlastpk1[i].bit = 6; goto SWITCH_OUT1; }
                sieve[idx0] |= BIT_MASK[6]; idx0 += delta0;
            case 7:
                if (idx0 >= maxi) { vlastpk1[i].bit = 7; goto SWITCH_OUT1; }
                sieve[idx0] |= BIT_MASK[7]; idx0 += delta1;
                };
        }
    SWITCH_OUT1: 
        vlastpk1[i].idx = idx0;
    }
    for (unsigned i = 0; i < maxi3; i++)
    {
        unsigned idx0 = vlastpk3[i].idx, delta0 = vlastpk3[i].delta, delta1 = delta0 + 1;
        switch (vlastpk3[i].bit)
        {
        default: NEVERHERE;
            for (;;)
            {
        case 0:
            if (idx0 >= maxi) { vlastpk3[i].bit = 0; goto SWITCH_OUT3; }
            sieve[idx0] |= BIT_MASK[0]; idx0 += delta0;
        case 3:
            if (idx0 >= maxi) { vlastpk3[i].bit = 3; goto SWITCH_OUT3; }
            sieve[idx0] |= BIT_MASK[3]; idx0 += delta0;
        case 6:
            if (idx0 >= maxi) { vlastpk3[i].bit = 6; goto SWITCH_OUT3; }
            sieve[idx0] |= BIT_MASK[6]; idx0 += delta1;
        case 1:
            if (idx0 >= maxi) { vlastpk3[i].bit = 1; goto SWITCH_OUT3; }
            sieve[idx0] |= BIT_MASK[1]; idx0 += delta0;
        case 4:
            if (idx0 >= maxi) { vlastpk3[i].bit = 4; goto SWITCH_OUT3; }
            sieve[idx0] |= BIT_MASK[4]; idx0 += delta0;
        case 7:
            if (idx0 >= maxi) { vlastpk3[i].bit = 7; goto SWITCH_OUT3; }
            sieve[idx0] |= BIT_MASK[7]; idx0 += delta1;
        case 2:
            if (idx0 >= maxi) { vlastpk3[i].bit = 2; goto SWITCH_OUT3; }
            sieve[idx0] |= BIT_MASK[2]; idx0 += delta0;
        case 5:
            if (idx0 >= maxi) { vlastpk3[i].bit = 5; goto SWITCH_OUT3; }
            sieve[idx0] |= BIT_MASK[5]; idx0 += delta1;
            }
        }
    SWITCH_OUT3:
        vlastpk3[i].idx = idx0;
    }
    for (unsigned i = 0; i < maxi5; i++)
    {
        unsigned idx0 = vlastpk5[i].idx, delta0 = vlastpk5[i].delta, delta1 = delta0 + 1;
        switch (vlastpk5[i].bit)
        {
        default: NEVERHERE;
            for (;;)
            {
        case 0:
            if (idx0 >= maxi) { vlastpk5[i].bit = 0; goto SWITCH_OUT5; }
            sieve[idx0] |= BIT_MASK[0]; idx0 += delta0;
        case 5:
            if (idx0 >= maxi) { vlastpk5[i].bit = 5; goto SWITCH_OUT5; }
            sieve[idx0] |= BIT_MASK[5]; idx0 += delta1;
        case 2:
            if (idx0 >= maxi) { vlastpk5[i].bit = 2; goto SWITCH_OUT5; }
            sieve[idx0] |= BIT_MASK[2]; idx0 += delta0;
        case 7:
            if (idx0 >= maxi) { vlastpk5[i].bit = 7; goto SWITCH_OUT5; }
            sieve[idx0] |= BIT_MASK[7]; idx0 += delta1;
        case 4:
            if (idx0 >= maxi) { vlastpk5[i].bit = 4; goto SWITCH_OUT5; }
            sieve[idx0] |= BIT_MASK[4]; idx0 += delta1;
        case 1:
            if (idx0 >= maxi) { vlastpk5[i].bit = 1; goto SWITCH_OUT5; }
            sieve[idx0] |= BIT_MASK[1]; idx0 += delta0;
        case 6:
            if (idx0 >= maxi) { vlastpk5[i].bit = 6; goto SWITCH_OUT5; }
            sieve[idx0] |= BIT_MASK[6]; idx0 += delta1;
        case 3:
            if (idx0 >= maxi) { vlastpk5[i].bit = 3; goto SWITCH_OUT5; }
            sieve[idx0] |= BIT_MASK[3]; idx0 += delta1;
            }
        }
    SWITCH_OUT5:
        vlastpk5[i].idx = idx0;
    }
    for (unsigned i = 0; i < maxi7; i++)
    {
        unsigned idx0 = vlastpk7[i].idx, delta0 = vlastpk7[i].delta, delta1 = delta0 + 1;
        switch (vlastpk7[i].bit)
        {
        default: NEVERHERE;
            for (;;)
            {
        case 7:
            if (idx0 >= maxi) { vlastpk7[i].bit = 7; goto SWITCH_OUT7; }
            sieve[idx0] |= BIT_MASK[7]; idx0 += delta1;
        case 6:
            if (idx0 >= maxi) { vlastpk7[i].bit = 6; goto SWITCH_OUT7; }
            sieve[idx0] |= BIT_MASK[6]; idx0 += delta1;
        case 5:
            if (idx0 >= maxi) { vlastpk7[i].bit = 5; goto SWITCH_OUT7; }
            sieve[idx0] |= BIT_MASK[5]; idx0 += delta1;
        case 4:
            if (idx0 >= maxi) { vlastpk7[i].bit = 4; goto SWITCH_OUT7; }
            sieve[idx0] |= BIT_MASK[4]; idx0 += delta1;
        case 3:
            if (idx0 >= maxi) { vlastpk7[i].bit = 3; goto SWITCH_OUT7; }
            sieve[idx0] |= BIT_MASK[3]; idx0 += delta1;
        case 2:
            if (idx0 >= maxi) { vlastpk7[i].bit = 2; goto SWITCH_OUT7; }
            sieve[idx0] |= BIT_MASK[2]; idx0 += delta1;
        case 1:
            if (idx0 >= maxi) { vlastpk7[i].bit = 1; goto SWITCH_OUT7; }
            sieve[idx0] |= BIT_MASK[1]; idx0 += delta1;
        case 0:
            if (idx0 >= maxi) { vlastpk7[i].bit = 0; goto SWITCH_OUT7; }
            sieve[idx0] |= BIT_MASK[0]; idx0 += delta0;
            }
        }
    SWITCH_OUT7:
        vlastpk7[i].idx = idx0;
    }
}

void InitLastPKB1(const uint64_t start, const uint64_t stop, const unsigned ridx)
{
    unsigned sml, med, lrg; sml = med = lrg = 0;
    maxi1 = maxi3 = maxi5 = maxi7 = 0;

    const uint64_t maxp = (uint64_t)floor(sqrt(stop));
    uint64_t p = first_root_prime;
    unsigned i = first_root_prime_idx;
    while (i < NO_ROOT_PRIMES)
    {   // for each root prime
        // get initial multiplier approximation
        const uint64_t startk = (start > p * p) ? (start / p) : p;
        // adjust to multiple of wheel_span
        uint64_t k = wheel.wheel_span * (startk / wheel.wheel_span);
        if (p < chunk_span)     sml++; else
        if (p < segment_span)   med++; else
                                lrg++;
        // move to current factor
        k += wheel.factor[rootsieveR[i]][ridx];
        if (k < startk) k += wheel.wheel_span;
        uint64_t n = k * p;
        // adjust if lower
        while (n < start) { k += wheel.wheel_span; n = k * p; };
        n = (n - start) / wheel.wheel_span;
        lastpk* pvals;
        switch (p % 8)
        {
        default: NEVERHERE;
        case 1: pvals = vlastpk1 + maxi1++; break;
        case 3: pvals = vlastpk3 + maxi3++; break;
        case 5: pvals = vlastpk5 + maxi5++; break;
        case 7: pvals = vlastpk7 + maxi7++; break;
        };
        //if (p % 8 > 3)
        //    if (p % 8 == 5)
        //        pvals = vlastpk5 + maxi5++;
        //    else   // == 7
        //        pvals = vlastpk7 + maxi7++;
        //else
        //    if (p % 8 == 1)
        //        pvals = vlastpk1 + maxi1++;
        //    else   // == 3
        //        pvals = vlastpk3 + maxi3++;
        assert(n / 8 < std::numeric_limits<unsigned>::max());
        pvals->idx = unsigned(n / 8); 
        pvals->bit = n % 8; 
        pvals->delta = unsigned(p / 8);

        // get next root prime value
        p += root_primes_gap[++i]; if (p > maxp) break;
    }
    //std::cout << maxi1 << " | " << maxi3 << " | " << maxi5 << " | " << maxi7; nln();

    assert(i < NO_ROOT_PRIMES); // should practically never reach the end of the list
    assert(med == 0 and lrg == 0);
}

uint64_t ProcessBatchB1(const uint64_t Nmax, const unsigned batch)
{
    const uint64_t segment = batch / wheel.wheel_values;
    const unsigned ridx = batch % wheel.wheel_values;
    const uint64_t seg_start = segment * segment_span; if (seg_start > Nmax) return 0;
    const uint64_t seg_stop = std::min(Nmax, seg_start + segment_span);
    const unsigned value = wheel.wheel[ridx];

    InitLastPKB1(seg_start, seg_stop, ridx);
    memset(sieve, 0, segment_bytes);
    //std::cout << imax << " | " << segment_bytes << " | " << segment_span; nln();

    for (uint64_t ck_start = seg_start; ck_start < seg_stop; ck_start += chunk_span)
    {
        const uint64_t ck_stop = std::min(ck_start + chunk_span, seg_stop);
        SieveChunkB1(seg_start, ck_stop, value);
    }

    uint64_t numPrimes = 0;
    const unsigned imax = unsigned((seg_stop - seg_start) / byte_bits / wheel.wheel_span);
    for (unsigned i = 0; i < imax; i++)
        numPrimes += bit0_count_table[sieve[i]];
    if (imax < segment_bytes)
        for (unsigned b = 0; b < byte_bits; b++)
            if (!(sieve[imax] & BIT_MASK[b]))
            {
                uint64_t n = seg_start + wheel.wheel_span * ((uint64_t)byte_bits * imax + b) + value;
                if (n > Nmax) break;
                numPrimes++;
            }
    //std::cout << seg_start << " | " << seg_stop << " | " << value; nln();
    return numPrimes;
}

uint64_t ProcessThreadB1(const uint64_t Nmax)
{
    // for each segment|value pair we have one batch 
    const uint64_t batches_no = segments_no * wheel.wheel_values;

    //vlastpk = new tpPrime[NO_ROOT_PRIMES]; vlastpk[0] = 0;
    vlastpk1 = new lastpk[szi1]; 
    vlastpk3 = new lastpk[szi3];
    vlastpk5 = new lastpk[szi5]; 
    vlastpk7 = new lastpk[szi7];
    sieve = new uint8_t[segment_bytes];
    auto _kk = __STDCPP_DEFAULT_NEW_ALIGNMENT__;
    //std::cout << std::hex << (uint64_t)(sieve) << " | "; nln();

    uint64_t numPrimes = 0;
    for (unsigned batch = NextBatchT(); batch < batches_no; batch = NextBatchT())
    { // do one batch
        numPrimes += ProcessBatchB1(Nmax, batch);
        fprintf(stdout, "%5.1f%%\b\b\b\b\b\b", NextBatchP() * 100.0 / batches_no);
    }
    delete[] vlastpk1; delete[] vlastpk3;
    delete[] vlastpk5; delete[] vlastpk7;
    delete[] sieve;

    return numPrimes;
}

uint64_t ParallelSieveB1(const uint64_t Nmax)
{
    const unsigned nativeThreads = std::max(1u, std::thread::hardware_concurrency());

    ResetSeed();
    SetSegmentsNo(Nmax, nativeThreads);

    uint64_t numPrimes = first_root_prime_idx + 1ull;
    std::vector<std::future<uint64_t>> counters;
    for (unsigned t = 0; t < nativeThreads; t++)
    {
        counters.push_back(std::async(&ProcessThreadB1, Nmax));
        //numPrimes += ProcessThreadB1(Nmax);
    }

    for (auto& c : counters)
        numPrimes += c.get();

    return numPrimes;
}

uint64_t SoE_bitpat_free1(const uint64_t Nmax, void*, void*, void*)
{
    InitRootPrimes();
    return ParallelSieveB1(Nmax);
}
