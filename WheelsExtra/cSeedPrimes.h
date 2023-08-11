#pragma once
#include "generators.h"
#include <mutex>

//template<tpPrime limit, unsigned size>
//constexpr auto GenSdPrmLst(void)
//{
//    //add 1 slot if required
//    constexpr tpPrime sieve_size = limit / 30 + ((limit % 30) > 0);
//    std::array <tpSieveElement, sieve_size> sieve;
//    for (tpPrime i = 0; i < sieve_size; i++) sieve[i] = 0xFF;
//    //mark 1 and numbers outside limit as composite
//    RESETbt(sieve, 1);
//    for (tpPrime i = 0; i < limit % 30; i++) RESETbt(sieve, limit + i);
//
//    tpPrime n = 1ull + root_sieve_steps[0];
//    unsigned i = 1;
//    unsigned ix = 1;
//    do
//    {
//        if (GETbt(sieve, n))
//        {
//            tpPrime n1 = n * n;
//            while (root_sieve_values[ix] != n % 30)
//                if (++ix == 8) ix = 0;
//            unsigned j = ix;
//            do
//            {
//                RESETbt(sieve, n1);
//                n1 += n * root_sieve_steps[j];
//                IncrementRootIdx(j);
//            } while (n1 < limit);
//        }
//        n += root_sieve_steps[i];
//        IncrementRootIdx(i);
//    } while (n < ConstSqrt(limit));
//
//    tpPrime numPrimes = 0;
//
//    std::array<tpRoot, size+1> primes{};
//
//    n = 0;
//    for (i = 0; i < sieve_size; i++)
//    {
//        auto ProcessBit = [&](const unsigned char msk, const unsigned char r)
//        {
//            if (sieve[i] & msk)
//            {
//                tpPrime nk = n + r;
//                if(nk < limit)
//                    primes[numPrimes++] = nk;
//            }
//        };
//
//        ProcessBit(1, 1); ProcessBit(2, 7);
//        ProcessBit(4, 11); ProcessBit(8, 13);
//        ProcessBit(16, 17); ProcessBit(32, 19);
//        ProcessBit(64, 23); ProcessBit(128, 29);
//
//        n += root_primorial;
//        if (numPrimes > 3300)
//            n = n;
//    };
//
//    return primes;
//}

constexpr unsigned LINES_PER_WINDOW = 9000;
constexpr unsigned BYTES_PER_WINDOW = LINES_PER_WINDOW * CACHE_LINE_BYTES;
constexpr unsigned CANDIDATES_PER_WINDOW = BYTES_PER_WINDOW * root_primorial;


struct cSeedPrimeElement
{
public:
    uint8_t prime_step = 0;
    uint8_t root_idx = 0;   // next root idx to be used
    tpPrime next_pos = 0;  // position in next window
#ifdef _DEBUG
    tpPrime prime_value = 0; 
#endif 

};

constexpr tpPrime prms_val_lmt[] = {
    ConstSqrt(1'000'000'000ull) ,
    ConstSqrt(10'000'000'000ull),
    ConstSqrt(50'000'000'000ull),
    ConstSqrt(100'000'000'000ull)
};
constexpr unsigned prms_no_lmt[] = {
    3'401 - root_N,
    9'592 - root_N,
    19'910 - root_N,
    27'293 - root_N
};

constexpr unsigned LMT_SMALL = 1u;
constexpr unsigned LMT_BIG = 3u;

#define WORK_BIG
#ifdef WORK_BIG
constexpr unsigned LMT_USED = LMT_BIG;
#else
constexpr unsigned LMT_USED = LMT_SMALL;
#endif

constexpr unsigned SEED_SIZE = prms_no_lmt[LMT_USED];
typedef std::array<cSeedPrimeElement, SEED_SIZE> cSeedElements;
typedef std::array<tpSieveElement, BYTES_PER_WINDOW> cSieveWindow;


class cSynchronizer
{
    std::mutex mtx;
    std::condition_variable cv;
    unsigned counter = 0;

public:
    inline void WaitTurnAndAcquire(unsigned turn)
    {
        std::unique_lock<std::mutex> lck(mtx);
        cv.wait(lck, [&] { return counter == turn; });
    }

    inline void Release(void)
    {
        {
            std::lock_guard<std::mutex> lk(mtx);
            counter++;
        };
        cv.notify_all();
    };
};

class cSeedPrimesGen
{
    CACHE_ALIGN cSeedElements seed_elements{};
    CACHE_ALIGN cSieveWindow sieve_window_unq{};
    
    CACHE_ALIGN cSieveWindow sieve_window_1{};
    CACHE_ALIGN cSieveWindow sieve_window_2{};
    cSynchronizer sync_1, sync_2;

    //std::mutex mtx_1, mtx_2;
    //std::condition_variable cv_1, cv_2;
    //unsigned window_state_1 = 0;
    //unsigned window_state_2 = 0;
    //std::atomic<unsigned> window_state_1 = 0;
    //std::atomic<unsigned> window_state_2 = 0;

    tpPrime limit = 0;
    tpPrime numPrimes = 0;

    cChecker checker;

    tpPrime CountWindowSimple(cSieveWindow& sieve_window, tpPrime window_start, unsigned win_sz);
    tpPrime CountWindow(cSieveWindow& sieve_window, tpPrime window_start);
    void MarkWindow(cSieveWindow& sieve_window, const tpPrime window_start);

    void InitializeSeedElements(void);
    void CountAll(void);
    void MarkAll(void);

public:
    tpPrime GenAllIn(const tpPrime lmt, bool dogen);
    tpPrime GenAllInPrl(const tpPrime lmt);
};
