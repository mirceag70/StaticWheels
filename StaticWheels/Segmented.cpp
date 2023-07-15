#include "Helper.h"

void SoA_LP_gen_root_primes(void);
extern uint8_t rsvR[NO_ROOT_PRIMES];
void GetRootPrimes(void);

consteval unsigned gcd(unsigned a, unsigned b) { return (a == 0) ? b : gcd(b % a, a); }
consteval unsigned lcm(unsigned a, unsigned b) { return (a * b) / gcd(a, b); }

constexpr uint16_t first_primes[] = 
        {2,         3,         5,         7,        11,        13,        17,        19,
		23,        29,        31,        37,        41,        43,        47,        53};
constexpr uint16_t wheel_entries[] = 
        { 0, 1, 2, 8, 48, 480, 5'760};

consteval unsigned Primorial(const unsigned n)
{
    unsigned prim = 1;
    for (unsigned i = 0; i < n; i++) prim *= first_primes[i];
    return prim;
}

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
static_assert(one_chunk_tail == 0); // only valid for wheel 3...
//constexpr unsigned min_chunk_span = lcm(8, wheel_values);
constexpr unsigned chunk_span = one_chunk_wheels * wheel_span;

constexpr unsigned chunk_size_bits = (1 << 14) * 8;

void SieveOneSeg(const uint64_t Nmax, const unsigned ridx, const unsigned no_t, const unsigned crnt_t)
{
    const tpPrime segsz = Nmax / no_t;
    const tpPrime start = crnt_t * segsz;
    const tpPrime stop = (crnt_t == (no_t - 1)) ? Nmax : start + segsz;

    //memset(sv, 0, 8 * svlen);
}

enum class defWorkType { parallel = true, sequential = false };

template<defWorkType workType>
void ParallelSegSieve(const uint64_t Nmax)
{
    static_assert(rootN == 3);  // only this one implemented here...

    constexpr unsigned nt = 3;
    constexpr std::array<const unsigned, wheel_values> vnt =
        { /*1*/nt + 0, /* 7*/nt + 0, /*11*/nt + 0, /*13*/nt + 0,
         /*17*/nt + 0, /*19*/nt + 0, /*23*/nt + 0, /*29*/nt + 0 };
    std::array<std::vector<std::thread>, wheel_values> vt;

    for (auto ridx : std::views::iota(0u, wheel_values))
        for (auto j : std::views::iota(0u, vnt[ridx]))
        {
            if constexpr (workType == defWorkType::parallel)
                vt[ridx].push_back(std::thread(&SieveOneSeg, Nmax, ridx, vnt[ridx], j));
            else
                SieveOneSeg( Nmax, ridx, vnt[ridx], j);
        }
    for (auto i : std::views::iota(0u, wheel_values)) for (auto& t : vt[i]) t.join();
}

uint64_t SoE_w8pat_seg(const uint64_t Nmax, void*, void*, void*)
{
    cTimer tmr; tmr.Start();
    GetRootPrimes();
    tmr.LapTime(true, "roots");

    ParallelSegSieve<defWorkType::parallel>(Nmax);

    //return CountAllSegments(vPrimes, Nmax, svlen);
    //return CountAllOnce(vPrimes, Nmax, svlen);

    tmr.Stop();
    return 0;
}
