#include "Helper.h"
#include "wheels.h"

constexpr unsigned WHEEL_N = 5; // [3..6]
wheels::RootSieve<WHEEL_N> wheel;

// bits in a byte, allways EIGHT
constexpr unsigned byte_bits = 8;
static_assert(byte_bits == 8); 

// the batch is divided in chunks
// chuncks in one batch are processed sequentially / incrementally
// dat for one chunk should fit in cache
constexpr unsigned chunk_bytes_lg2 = 16;   // 2^15 = 32k
constexpr unsigned chunk_bytes = 1 << chunk_bytes_lg2;  
// try to have sizes as 2^n everytime
static_assert(chunk_bytes == 1 << chunk_bytes_lg2);
constexpr uint64_t chunk_bits = chunk_bytes * byte_bits;
constexpr uint64_t chunk_span = chunk_bits * wheel.wheel_span;

// the whole interval is divided in segments
// each segment should fit N full chunks, with no rest
//constexpr unsigned max_segment_chunks_lg2 = 13;   // 2^13 = 8k
//constexpr unsigned max_segment_chunks = 1 << max_segment_chunks_lg2;
//static_assert(max_segment_chunks == (1 << max_segment_chunks_lg2));
//static_assert(max_segment_chunks * chunk_bytes < std::numeric_limits<unsigned>::max());
//constexpr unsigned max_segment_bytes = max_segment_chunks * chunk_bytes;
//static_assert(max_segment_bytes < std::numeric_limits<std::size_t>::max()/*SIZE_MAX*/);
constexpr unsigned max_segment_bytes = 2 << 30; //2GB


void SoA_LP_gen_root_primes(void);
extern uint8_t root_primes_gap[];
extern uint8_t rsvR[];
typedef wheels::RootSieve<WHEEL_N>::tpWheelEntry rootsieveRtp;
rootsieveRtp rootsieveR[NO_ROOT_PRIMES];
unsigned first_root_prime = wheel.wheel[1];
unsigned first_root_prime_idx;
void InitRootPrimes(void)
{
    SoA_LP_gen_root_primes();
    unsigned p = 5, i = 1;
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
    //const uint64_t maxp = (uint64_t)floor(sqrt(stop));
    uint64_t p = first_root_prime;
    unsigned i = first_root_prime_idx;
    while (i < NO_ROOT_PRIMES)
    {   // for each root prime
        // get initial multiplier approximation
        const uint64_t startk = (start > p * p) ? (start / p) : p;
        // adjust to multiple of wheel_span
        uint64_t k = wheel.wheel_span * (startk / wheel.wheel_span);
        // move to current factor
        k += wheel.factor[rootsieveR[i]][ridx];
        if (k < startk) k += wheel.wheel_span;
        uint64_t n = k * p;
        // adjust if lower
        while (n < start) { k += wheel.wheel_span; n = k * p; }
        vlastpk[i] = (n - start) / wheel.wheel_span;
        // get next root prime value
        p += root_primes_gap[++i]; if (p * p > stop) break;
    }
    assert(i < NO_ROOT_PRIMES); // should practically never reach the end of the list
}

uint64_t segments_no, segment_bytes, segment_span;

void SieveChunk(const uint64_t seg_start, const uint64_t stop, const uint64_t val)
{
    //const uint64_t maxn = (stop - start) / wheel.wheel_span;
    const uint64_t maxn = std::min (segment_bytes * byte_bits, (stop - seg_start) / wheel.wheel_span + 1);
    const uint64_t maxp = (uint64_t)floor(sqrt(stop));
    uint64_t p = first_root_prime;
    for (unsigned i = first_root_prime_idx; i < NO_ROOT_PRIMES;)
    {
        tpPrime n = vlastpk[i];
        if (n < maxn)
        {
            //uint64_t idx0 = n / 8, bit0 = n % 8;
            //uint64_t idxdeltap = p / 8, idxdeltab = p % 8;
            for (; n < maxn; n += p)
            {
                uint64_t idx = n / 8, bit = n % 8;
                uint64_t nn = seg_start + wheel.wheel_span * (8 * idx + bit) + val;
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
        p += root_primes_gap[++i]; if (p > maxp) break;
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
    //memset(sieve, 0, imax + 1ull);
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
    const uint64_t batches_no = segments_no * wheel.wheel_values;

    vlastpk = new tpPrime[NO_ROOT_PRIMES]; vlastpk[0] = 0;
    sieve = new uint8_t[segment_bytes];
    //std::cout << std::hex << (uint64_t)(sieve) << " | "; nln();

    uint64_t numPrimes = 0;
    for (unsigned batch = NextBatchT(); batch < batches_no; batch = NextBatchT())
    { // do one batch
        numPrimes += ProcessBatch(Nmax, batch);
        fprintf(stdout, "%4.1f%%\b\b\b\b\b", NextBatchP() * 100.0 / batches_no);
    }
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
    std::cout << " segments; " << wheel.wheel_values << " values ]";
}

uint64_t ParallelSieve(const uint64_t Nmax)
{
    const unsigned nativeThreads = std::max(1u, std::thread::hardware_concurrency());

    ResetSeed();
    SetSegmentsNo(Nmax, nativeThreads);

    uint64_t numPrimes = first_root_prime_idx + 1ull;
    std::vector<std::future<uint64_t>> counters;
    for (unsigned t = 0; t < nativeThreads; t++)
    {
        counters.push_back(std::async(&ProcessThread, Nmax));
        //numPrimes += ProcessThread(Nmax);
    }
    fprintf(stdout, "     ");

    for (auto& c : counters) 
        numPrimes += c.get();

    return numPrimes;
}

uint64_t SoE_pat_free(const uint64_t Nmax, void*, void*, void*)
{
    InitRootPrimes();
    return ParallelSieve(Nmax);
}
