#pragma once
#include "Helper.h"
#include <mutex>

constexpr auto CACHE_LINE_BYTES = 32;
#define CACHE_ALIGN __declspec(align(CACHE_LINE_BYTES))

typedef  uint8_t tpSieveElement;
static_assert(sizeof(tpSieveElement) == 1);
template <size_t SZ>
using tpSieveArray = std::array<tpSieveElement, SZ>;

consteval uint8_t Msk(unsigned i) { return (1 << (i)); };

constexpr unsigned root_N = 3;
constexpr unsigned root_sieve_size = 8;
constexpr unsigned root_primorial = 2 * 3 * 5;
constexpr CACHE_ALIGN tpSieveArray<root_sieve_size> root_sieve_values
        { 1, 7, 11, 13, 17, 19, 23, 29 };
constexpr CACHE_ALIGN tpSieveArray<root_sieve_size> root_sieve_steps
		{ 6, 4, 2, 4, 2, 4, 6, 2 };
constexpr CACHE_ALIGN tpSieveArray<root_sieve_size> bit_masks
        { Msk(0), Msk(1), Msk(2), Msk(3), Msk(4), Msk(5), Msk(6), Msk(7) };
constexpr CACHE_ALIGN tpSieveArray<root_sieve_size> bit_inverse_masks
        { ~Msk(0), ~Msk(1), ~Msk(2), ~Msk(3), ~Msk(4), ~Msk(5), ~Msk(6), ~Msk(7) };
constexpr CACHE_ALIGN tpSieveArray<root_primorial> root_sieve_masks
        //{ 0, Msk(0), 0, 0, 0,
        //  0, 0, Msk(1), 0, 0,
        //  0, Msk(2), 0, Msk(3), 0,
        //  0, 0, Msk(4), 0, Msk(5),
        //  0, 0, 0, Msk(6), 0,
        //  0, 0, 0, 0, Msk(7) };
        {200, Msk(0), 202, 203, 204,
         205, 206, Msk(1), 208, 209,
         210, Msk(2), 212, Msk(3), 214,
         215, 216, Msk(4), 218, Msk(5),
         220, 221, 222, Msk(6), 224,
         225, 226, 227, 228, Msk(7) };

template<typename T1, typename T2>
constexpr bool GETbt(T1& vect, const T2 pos)
{
    T2 q = pos / 30;
    T2 r = pos % 30;
    tpSieveElement msk = root_sieve_masks[r];
    //assert(msk < 200);
    return (vect[q] & msk);
}

template<typename T1, typename T2>
constexpr void RESETbt(T1& vect, const T2 pos)
{
    T2 q = pos / 30;
    T2 r = pos % 30;
    tpSieveElement msk = root_sieve_masks[r];
    //assert(msk < 200);
    vect[q] &= ~msk;
}

constexpr void IncrementRootIdx(uint8_t& idx)
{
    idx == (root_sieve_size - 1) ? idx = 0 : idx++;
}

tpRoot TestPrimes(void);

constexpr unsigned LINES_PER_WINDOW = 9000;
constexpr unsigned BYTES_PER_WINDOW = LINES_PER_WINDOW * CACHE_LINE_BYTES;
constexpr unsigned BITS_PER_WINDOW = BYTES_PER_WINDOW * 8;
constexpr unsigned WINDOW_SIZE = BITS_PER_WINDOW * root_primorial;

constexpr unsigned LMT_SMALL = 1u;
constexpr unsigned LMT_BIG = 3u;

#define WORK_BIG
#ifdef WORK_BIG
constexpr unsigned LMT_USED = LMT_BIG;
#else
constexpr unsigned LMT_USED = LMT_SMALL;
#endif

constexpr CACHE_ALIGN tpPrime prms_val_lmt[]
        {
            ConstSqrt(1'000'000'000ull) ,
            ConstSqrt(10'000'000'000ull),
            ConstSqrt(50'000'000'000ull),
            ConstSqrt(100'000'000'000ull),
            ConstSqrt(1'000'000'000'000ull)
        };
constexpr CACHE_ALIGN unsigned prms_no_lmt[]
        {
            3'401 - root_N,
            9'592 - root_N,
            19'910 - root_N,
            27'293 - root_N,
            270'293 - root_N
        };

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

struct cSeedPrimeElement
{
public:
    uint8_t prime_step = 0;
    tpPrime next_pos = 0;  // position in next window
#ifdef _DEBUG
    tpPrime prime_value = 0;
#endif 
};

constexpr unsigned SEED_SIZE = prms_no_lmt[LMT_USED];
typedef std::array<cSeedPrimeElement, SEED_SIZE> cSeedElements;
typedef tpSieveArray<BYTES_PER_WINDOW> cSieveWindow;

struct cAllWindows
{
    cSieveWindow CACHE_ALIGN window_01{};
    cSieveWindow CACHE_ALIGN window_07{};
    cSieveWindow CACHE_ALIGN window_11{};
    cSieveWindow CACHE_ALIGN window_13{};
    cSieveWindow CACHE_ALIGN window_17{};
    cSieveWindow CACHE_ALIGN window_19{};
    cSieveWindow CACHE_ALIGN window_23{};
    cSieveWindow CACHE_ALIGN window_29{};
};

struct cAllSeedElements
{
    cSeedElements CACHE_ALIGN seeds_01{};
    cSeedElements CACHE_ALIGN seeds_07{};
    cSeedElements CACHE_ALIGN seeds_11{};
    cSeedElements CACHE_ALIGN seeds_13{};
    cSeedElements CACHE_ALIGN seeds_17{};
    cSeedElements CACHE_ALIGN seeds_19{};
    cSeedElements CACHE_ALIGN seeds_23{};
    cSeedElements CACHE_ALIGN seeds_29{};
};

class cFastGen
{
	tpPrime limit = 0;
    tpPrime num_primes = 0;

    cChecker checker;

    static inline cAllSeedElements all_seed_elements{};
    static inline cAllWindows all_sieve_windows1{};
    static inline cAllWindows all_sieve_windows2{};
    cSynchronizer sync_1, sync_2;

    typedef std::array<cSieveWindow*, 8> svsarr;

    static inline svsarr wndws1 = 
      { & all_sieve_windows1.window_01, & all_sieve_windows1.window_07,
        & all_sieve_windows1.window_11, & all_sieve_windows1.window_13, 
        & all_sieve_windows1.window_17, & all_sieve_windows1.window_19, 
        & all_sieve_windows1.window_23, & all_sieve_windows1.window_29 };

    static inline svsarr wndws2 =
    { &all_sieve_windows2.window_01, &all_sieve_windows2.window_07,
      &all_sieve_windows2.window_11, &all_sieve_windows2.window_13,
      &all_sieve_windows2.window_17, &all_sieve_windows2.window_19,
      &all_sieve_windows2.window_23, &all_sieve_windows2.window_29 };


    void InitializeSeedElements(cSeedElements& seed_elements, const unsigned rem);

    void MarkWindow(cSieveWindow& window, cSeedElements& elements);
    tpPrime CountWindows(cAllWindows& windows, const tpPrime window_start);
    tpPrime CountWindowsSimple(cAllWindows& windows, const tpPrime window_start, const unsigned win_sz);
    tpPrime CountWindow(cSieveWindow& sieve_window, const uint8_t r, const tpPrime window_start);


    void MarkTT(void);
    void CountTT(void);

    void MarkNCount(void);

    tpPrime CountWindowsIterate(cAllWindows& windows, const tpPrime window_start);

public:
    cFastGen(void);
   
    tpPrime GenPrimes(const tpPrime lmt, bool dogen);
    tpPrime GenAllInPrl(const tpPrime lmt);
    //tpPrime GenAllNotPrl(const tpPrime lmt);
};

