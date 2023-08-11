#include "cSeedPrimes.h"


template<unsigned LMT>
#ifdef WORK_BIG
constexpr auto GenerateSeedPrimesList(void)
#else
consteval auto GenerateSeedPrimesList(void)
#endif
{
    static_assert(LMT >= 0 and LMT < 4);
    constexpr tpPrime limit = prms_val_lmt[LMT];
    constexpr unsigned tail = limit % root_primorial;
    //add 1 slot if required
    constexpr tpPrime sieve_size = limit / root_primorial + (tail > 0);
    std::array<tpSieveElement, sieve_size> sieve; sieve.fill(0xFF);
    //mark 1 and numbers outside limit as composite
    RESETbt(sieve, 1);
    if (tail > 0)
        for (uint8_t i = 0; i < (root_primorial - tail); i++)
        {
            switch ((limit + i) % 30)
            {
            case  1: case  7: case 11: case 13:
            case 17: case 19: case 23: case 29:
                RESETbt(sieve, limit + i);
            }
        }

    tpPrime n = 1ull + root_sieve_steps[0];
    uint8_t i = 1;
    uint8_t ix = 1;
    do
    {
        if (GETbt(sieve, n))
        {
            tpPrime n1 = n * n;
            while (root_sieve_values[ix] != n % root_primorial)
                IncrementRootIdx(ix);
            uint8_t j = ix;
            do
            {
                RESETbt(sieve, n1);
#pragma warning(suppress:28020) // expression not true at this call
                n1 += n * root_sieve_steps[j];
                IncrementRootIdx(j);
            } while (n1 < limit);
        }
        n += root_sieve_steps[i];
        IncrementRootIdx(i);
    } 
    while (n < ConstSqrt(limit));

    constexpr unsigned size = prms_no_lmt[LMT];
    std::array<tpRoot, size> primes{};

    tpPrime numPrimes = 0;
    unsigned k;
    for (n = k = 0; k < sieve_size; k++)
    {
        auto ProcessBit = [&](const unsigned char msk, const unsigned char r)
        {
            if (sieve[k] & msk)
            {
                primes[numPrimes++] = (tpRoot)(n + r);
            }
        };

        ProcessBit(1, 1); ProcessBit(2, 7);
        ProcessBit(4, 11); ProcessBit(8, 13);
        ProcessBit(16, 17); ProcessBit(32, 19);
        ProcessBit(64, 23); ProcessBit(128, 29);

        n += root_primorial;
    };

    return primes;
}

#ifdef WORK_BIG
const auto CACHE_ALIGN seed_primes = GenerateSeedPrimesList<LMT_USED>();
#else
constexpr auto CACHE_ALIGN seed_primes = GenerateSeedPrimesList<LMT_USED>();
#endif

//tpRoot TestPrimes(void)
//{
//    cChecker checker;
//    unsigned i;
//    for (i = 0; i < root_N; i++)
//    {
//        assert(checker.check_next_prime(firstPrimes[i]));
//    }
//    for (auto n : seed_primes)
//    {
//        assert(checker.check_next_prime(n));
//        i++;
//    }
//    return 0;
//}

//constexpr uint8_t GetBit(unsigned rem)
//{
//    switch (rem)
//    {
//    case  1: return 0;
//    case  7: return 1;
//    case 11: return 2;
//    case 13: return 3;
//    case 17: return 4;
//    case 19: return 5;
//    case 23: return 6;
//    case 29: return 7;
//    default:
//        assert(false);
//        return 0;
//    }
//}
//
//constexpr void ComputeNextWindow(cSeedPrimeElement& element, tpPrime val)
//{
//    assert((val / CANDIDATES_PER_WINDOW) < 0xFFFF'FFFF);
//    element.next_win = (uint32_t)(val / CANDIDATES_PER_WINDOW);
//    val = val % CANDIDATES_PER_WINDOW;
//    assert((val / root_primorial) < 0xFFFF);
//    element.next_win_pos = (uint16_t)val;
//}

template<unsigned LMT>
consteval auto GenerateSeedElements(void)
//constexpr auto GenerateSeedElements(void)
{
    constexpr unsigned size = prms_no_lmt[LMT];
    std::array<cSeedPrimeElement, size> elements{};
    uint8_t ix = 1;

    auto InitializeValues = [&](unsigned i)
    {
        //start from n*n
        elements[i].next_pos = seed_primes[i] * seed_primes[i];
#pragma warning(suppress:28020) // expression not true at this call
        while (root_sieve_values[ix] != (seed_primes[i] % root_primorial))
            IncrementRootIdx(ix);
        elements[i].root_idx = ix;
    };

    unsigned i;
    for (i = 0; i < size - 1; i++)
    {
        elements[i].prime_step = seed_primes[i + 1] - seed_primes[i];
        InitializeValues(i);
    }
#pragma warning(suppress:28020) // expression not true at this call
    elements[i].prime_step = root_primorial + seed_primes[0] - seed_primes[i];
    InitializeValues(i);
    return elements;
}

//auto CACHE_ALIGN seed_elements = GenerateSeedElements<LMT_USED>();

//__declspec(noinline)
void cSeedPrimesGen::InitializeSeedElements(void)
{
    uint8_t ix = 1;

    auto InitializeValues = [&](unsigned i)
    {
        //start from n*n
        seed_elements[i].next_pos = seed_primes[i];
        seed_elements[i].next_pos *= seed_primes[i];
#pragma warning(suppress:28020) // expression not true at this call
        while (root_sieve_values[ix] != (seed_primes[i] % root_primorial))
            IncrementRootIdx(ix);
        seed_elements[i].root_idx = ix;
    };

    unsigned i;
    for (i = 0; i < seed_elements.size() - 1; i++)
    {
#ifdef _DEBUG
        seed_elements[i].prime_value = seed_primes[i];
#endif
        seed_elements[i].prime_step = seed_primes[i + 1] - seed_primes[i];
        InitializeValues(i);
    }
#ifdef _DEBUG
    seed_elements[i].prime_value = seed_primes[i];
#endif
    seed_elements[i].prime_step = root_primorial + seed_primes[0] - seed_primes[i];
    InitializeValues(i);

    for (i = 0; i < root_N; i++)
        assert(checker.check_next_prime(firstPrimes[i]));
}

//__declspec(noinline)
//__forceinline
void cSeedPrimesGen::MarkWindow(cSieveWindow& sieve_window, tpPrime window_start)
{
    sieve_window.fill(0xFF);
    tpPrime n = firstPrimes[root_N];
    for (auto& element : seed_elements)
    {
        assert(n == element.prime_value);
        //tpPrime n = element.prime_value;
        tpPrime position = element.next_pos;
        uint8_t idx = element.root_idx;
        while (position < CANDIDATES_PER_WINDOW)
        {
            //mark composite
            RESETbt(sieve_window, position);
            //cntr++;
            //advance root
#pragma warning(suppress:28020) // expression not true at this call
            position += n * root_sieve_steps[idx];
            IncrementRootIdx(idx);
        }
        element.root_idx = idx;
        element.next_pos = position - CANDIDATES_PER_WINDOW;
        n += element.prime_step;
    }
}

//__forceinline
tpPrime cSeedPrimesGen::CountWindow(cSieveWindow& sieve_window, const tpPrime window_start)
{
    tpPrime numPrimes = 0;

    tpPrime n = window_start;
    tpPrime k = 0;

    auto ProcessBit = [&](const unsigned char msk, const unsigned char r)
    {
        if (sieve_window[k] & msk)
        {
            tpPrime nk = n + r;
            if (/*nk > 1 and */nk < limit)
            {

                assert(nk == 1 or checker.check_next_prime(nk));
                numPrimes++;
            }
        }
    };

    for (; k < BYTES_PER_WINDOW; k++)
    {
        ProcessBit(1, 1); ProcessBit(2, 7);
        ProcessBit(4, 11); ProcessBit(8, 13);
        ProcessBit(16, 17); ProcessBit(32, 19);
        ProcessBit(64, 23); ProcessBit(128, 29);

        n += 30;
    };

    return numPrimes;
}

tpPrime cSeedPrimesGen::CountWindowSimple(cSieveWindow& sieve_window, const tpPrime window_start, const unsigned win_sz)
{
    tpPrime numPrimes = 0, k = 0;
    for (; k < win_sz; k++)
    {
        numPrimes += bit1_count_table[sieve_window[k]];
    };

    tpPrime n;
    auto ProcessBit = [&](const unsigned char msk, const unsigned char r)
    {
        if (sieve_window[k] & msk)
        {
            tpPrime nk = n + r;
            if (nk < limit)
            {
                numPrimes++;
            }
        }
    };
    if (win_sz < BYTES_PER_WINDOW)
    {
        n = window_start + 30 * k;
        ProcessBit(1, 1); ProcessBit(2, 7);
        ProcessBit(4, 11); ProcessBit(8, 13);
        ProcessBit(16, 17); ProcessBit(32, 19);
        ProcessBit(64, 23); ProcessBit(128, 29);
    }

    return numPrimes;
}

//__declspec(noinline)
tpPrime cSeedPrimesGen::GenAllIn(const tpPrime lmt, bool dogen)
{
    limit = lmt;
    
    // 2, 3, 5 are not generated but 1 is
    tpPrime numPrimes = 3 - 1;

    //cTimer tmr; tmr.Start();
    InitializeSeedElements();
    //tmr.Stop(true);

    for(tpPrime window_start = 0; window_start < limit; 
        window_start += CANDIDATES_PER_WINDOW)
    {
        MarkWindow(sieve_window_unq, window_start);
        numPrimes += dogen ?
            CountWindow(sieve_window_unq, window_start) :
            CountWindowSimple(sieve_window_unq, window_start, 
                            std::min(BYTES_PER_WINDOW, unsigned((limit - window_start)/30)));
    }

    return numPrimes;
}

//__declspec(noinline)
tpPrime cSeedPrimesGen::GenAllInPrl(const tpPrime lmt)
{
    limit = lmt;

    // 2, 3 and 5 are not generated
    numPrimes = 3;

    InitializeSeedElements();

    std::thread tm(&cSeedPrimesGen::MarkAll, this);
    std::thread tc(&cSeedPrimesGen::CountAll, this);
    tm.join(); 
    tc.join();

    return numPrimes;
}

void cSeedPrimesGen::MarkAll(void)
{
    bool tictoc = true;
    tpPrime window_start_1 = 0;
    tpPrime window_start_2 = CANDIDATES_PER_WINDOW;
    unsigned expected_winstate_1 = 0;
    unsigned expected_winstate_2 = 0;

    do
    {
        if (tictoc)
        {
            if (window_start_1 >= limit)
                break;

            tictoc = false;

            sync_1.WaitTurnAndAcquire(expected_winstate_1);
            expected_winstate_1 += 2;

            MarkWindow(sieve_window_1, window_start_1);
            window_start_1 += 2 * CANDIDATES_PER_WINDOW;

            sync_1.Release();
        }
        else
        {
            if (window_start_2 >= limit)
                break;

            tictoc = true;

            sync_2.WaitTurnAndAcquire(expected_winstate_2);
            expected_winstate_2 += 2;

            MarkWindow(sieve_window_2, window_start_2);
            window_start_2 += 2 * CANDIDATES_PER_WINDOW;

            sync_2.Release();
        }
    } while (true); 
}

void cSeedPrimesGen::CountAll(void)
{
    bool tictoc = true;
    tpPrime window_start_1 = 0;
    tpPrime window_start_2 = CANDIDATES_PER_WINDOW;
    unsigned expected_winstate_1 = 1;
    unsigned expected_winstate_2 = 1;
    do
    {
        if (tictoc)
        {
            if (window_start_1 >= limit)
                break;

            tictoc = false;

            sync_1.WaitTurnAndAcquire(expected_winstate_1);
            expected_winstate_1 += 2;

            numPrimes += CountWindow(sieve_window_1, window_start_1);
            window_start_1 += 2 * CANDIDATES_PER_WINDOW;

            sync_1.Release();
        }
        else
        {
            if (window_start_2 >= limit)
                break;

            tictoc = true;

            sync_2.WaitTurnAndAcquire(expected_winstate_2);
            expected_winstate_2 += 2;

            numPrimes += CountWindow(sieve_window_2, window_start_2);
            window_start_2 += 2 * CANDIDATES_PER_WINDOW;

            sync_2.Release();
        }
    } while (true); 
}

void testCachedSieve(const uint64_t LIMIT, bool dogen)
{
    std::cout << "\n - Cached | Basic " 
        << (dogen ? "w.gen.primes" : "count only") << " - \n";

    cTimer tmr;
    std::vector<decltype(tmr.LapTime())> times;

    tmr.Start();

    for (auto i : Range<1>())
    {
        cSeedPrimesGen *worker = new cSeedPrimesGen;
        tpPrime n = worker->GenAllIn(LIMIT, dogen);  nln();
        std::cout << n << " primes up to " << LIMIT; nln();
        times.push_back(tmr.LapTime(true));
        
        delete worker;
    }

    tmr.Stop();
    nln(true);
    std::cout << "Average compute time: " << Average(times);
}

void testCachedSieve2T(const uint64_t LIMIT)
{
    std::cout << "\n - Cached | Basic - 2 threads - \n";

    cTimer tmr;
    std::vector<decltype(tmr.LapTime())> times;

    tmr.Start();

    for (auto i : Range<5>())
    {
        cSeedPrimesGen* worker = new cSeedPrimesGen;
        tpPrime n = worker->GenAllInPrl(LIMIT);  nln();
        std::cout << n << " primes up to " << LIMIT; nln();
        times.push_back(tmr.LapTime(true));

        delete worker;
    }

    tmr.Stop();
    nln(true);
    std::cout << "Average compute time: " << Average(times);
}