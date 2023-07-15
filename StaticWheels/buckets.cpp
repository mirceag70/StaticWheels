#include "Helper.h"


constexpr uint16_t first_primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53 };

consteval unsigned Primorial(const unsigned n)
{
    unsigned prim = 1;
    for (unsigned i = 0; i < n; i++) prim *= first_primes[i];
    return prim;
}

constexpr uint16_t wheel_entries[] = { 0, 1, 2, 8, 48, 480, 5'760 };


extern uint8_t root_primes_gap[NO_ROOT_PRIMES];
//void SoA_LP_gen_root_primes(void);
extern uint8_t rsvR[NO_ROOT_PRIMES];
void GetRootPrimes(void);


// used wheel 
constexpr unsigned rootN = 3; static_assert(rootN >= 3 and rootN <= 6);


// wheel length
constexpr unsigned wheel_span = Primorial(rootN);
constexpr unsigned wheel_values = wheel_entries[rootN];
// chunk sizes - cache friendly
constexpr unsigned chunk_bytes_lg2 = 14;
constexpr unsigned chunk_bytes = (1 << chunk_bytes_lg2);
constexpr unsigned chunk_bits = 8 * chunk_bytes;
constexpr unsigned one_chunk_wheels = chunk_bits / wheel_values;
constexpr unsigned one_chunk_tail = chunk_bits % wheel_values;
//static_assert(one_chunk_tail == 0); // only valid for wheel 3...
constexpr unsigned chunk_span = one_chunk_wheels * wheel_span;


constexpr uint8_t wheel_factors[wheel_values][wheel_values] = {};


void SieveOneChunk(const unsigned ridx, const uint64_t start, const uint64_t stop, uint8_t sv[], tpPrime lastnsp[])
{
    const uint64_t maxn = stop / 30;
    uint64_t last_p = 5;
    for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
    {
        uint64_t p = last_p + root_primes_gap[i]; if (p * p > stop) break;
        tpPrime n = lastnsp[i];
        if (n <= maxn)
        {
            for (; n <= maxn; n += p)
            {
                SetBit(n, sv);
            }
            lastnsp[i] = n;
        }
        last_p = p;
    }
}

void SieveOneSeg(const uint64_t Nmax, uint8_t sv[], const unsigned ridx, const unsigned no_t, const unsigned crnt_t)
{
    const tpPrime segsz = Nmax / no_t;
    const tpPrime start = crnt_t * segsz;
    const tpPrime stop = (crnt_t == (no_t - 1)) ? Nmax : start + segsz;

    //init multiples vector
    tpPrime* lastnsp = new tpPrime[NO_ROOT_PRIMES]; lastnsp[0] = 0;
    uint64_t p = 5;
    for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
    {
        p += root_primes_gap[i]; if (p * p > stop) break;
        const unsigned deltak = wheel_factors[rsvR[i]][ridx];
        tpPrime k = deltak;
        if (start == 0)
        {
            k += 30ull * (p / 30); if (k < p) k += 30;
        }
        else
        {
            k += 30ull * (start / p / 30); if (k == 1) k = 31;
        }
        const tpPrime n = k * p; lastnsp[i] = n / 30;
    }

    tpPrime iter_min = start;
    if (stop > start + chunk_bits)
        for (; iter_min < (stop - chunk_bits); iter_min += chunk_bits)
        {
            SieveOneChunk(ridx, iter_min, iter_min + chunk_bits, sv, lastnsp);
        }
    SieveOneChunk(ridx, iter_min, stop, sv, lastnsp);

    delete[] lastnsp;
}

enum class defWorkType { parallel = true, sequential = false };
template<defWorkType workType>
void SieveOneR(const uint64_t Nmax, uint8_t sv[], const unsigned ridx)
{
    static_assert(rootN == 3);  // only this wheel implemented here...

    const tpPrime rsvlen = Nmax / (wheel_span * 8ull) + 1;
    uint8_t *rsv = sv + ridx * rsvlen;
    memset((void*)rsv, 0, rsvlen);

    constexpr unsigned nt = 3;
    //constexpr std::array<const unsigned, wheel_values> vnt =
    //{ /*1*/nt + 0, /* 7*/nt + 0, /*11*/nt + 0, /*13*/nt + 0,
    // /*17*/nt + 0, /*19*/nt + 0, /*23*/nt + 0, /*29*/nt + 0 };
    std::vector<std::thread> vt;
    auto Sieve = [&](const unsigned i)
    {
        if constexpr (workType == defWorkType::parallel)
            vt.push_back(std::thread(&SieveOneSeg, Nmax, rsv, ridx, nt, i));
        else
            SieveOneSeg(Nmax, rsv, ridx, nt, i);
    };
    for (auto i : std::views::iota(0u, nt))
        Sieve(i);
    for (auto& t : vt) t.join();
}

template<defWorkType workType>
void ParallelSegSieve(const uint64_t Nmax, uint8_t sv[])
{
    std::vector<std::thread> vt;

    auto Sieve = [&](const unsigned ridx)
    {
        if constexpr (workType == defWorkType::parallel)
            vt.push_back(std::thread(&SieveOneR<workType>, Nmax, sv, ridx));
        else
            SieveOneR<workType>(Nmax, sv, ridx);
    };

    for (auto ridx : std::views::iota(0u, wheel_values)) 
        Sieve(ridx);
    for (auto& t : vt) t.join();
}

tpPrime CountSeg(const uint64_t Nmax, uint8_t sv[]) { return 55; }

uint64_t SoE_w8pat_bckt(const uint64_t Nmax, uint8_t sv[], void*, void*)
{
    cTimer tmr; tmr.Start();
    GetRootPrimes();
    tmr.LapTime(true, "roots");

    //ParallelSegSieve<defWorkType::parallel>(Nmax, sv);
    ParallelSegSieve<defWorkType::sequential>(Nmax, sv);
    tmr.LapTime(true, "sieve");

    tpPrime numprimes = CountSeg(Nmax, sv);
    //tpPrime numprimes = CountAllOnce(vPrimes, Nmax, svlen);
    tmr.LapTime(true, "count");

    tmr.Stop();
    return numprimes;
}

