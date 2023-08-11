#include "cFastGen.h"

template<unsigned LMT>
#ifdef WORK_BIG
constexpr auto GenerateSeedPrimesList(void)
#else
consteval auto GenerateSeedPrimesList(void)
#endif
{
    static_assert(LMT >= 0 and LMT <= 3);
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
                n1 += n * root_sieve_steps[j];
                IncrementRootIdx(j);
            } while (n1 < limit);
        }
        n += root_sieve_steps[i];
        IncrementRootIdx(i);
    } while (n < ConstSqrt(limit));

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

cFastGen::cFastGen(void)
{
    InitializeSeedElements(all_seed_elements.seeds_01, 1);
    InitializeSeedElements(all_seed_elements.seeds_07, 7);
    InitializeSeedElements(all_seed_elements.seeds_11, 11);
    InitializeSeedElements(all_seed_elements.seeds_13, 13);
    InitializeSeedElements(all_seed_elements.seeds_17, 17);
    InitializeSeedElements(all_seed_elements.seeds_19, 19);
    InitializeSeedElements(all_seed_elements.seeds_23, 23);
    InitializeSeedElements(all_seed_elements.seeds_29, 29);

    for (int i = 0; i < root_N; i++)
        assert(checker.check_next_prime(firstPrimes[i]));
}

//__declspec(noinline)
void cFastGen::InitializeSeedElements(cSeedElements& seed_elements, const unsigned rem)
{
    auto InitializeValues = [&](int i)
    {
#ifdef _DEBUG
        seed_elements[i].prime_value = seed_primes[i];
#endif
        //start from n*n
        tpPrime first_val = seed_primes[i];
        first_val *= seed_primes[i];
        //advance to the first pos for rem
        while ((first_val % root_primorial) != rem)
            first_val += seed_primes[i];
        seed_elements[i].next_pos = first_val / root_primorial;
    };

    int i;
    for (i = 0; i < seed_elements.size() - 1; i++)
    {
        seed_elements[i].prime_step = seed_primes[i + 1] - seed_primes[i];
        InitializeValues(i);
    }
#pragma warning(suppress:28020) // expression not true at this call
    seed_elements[i].prime_step = root_primorial + seed_primes[0] - seed_primes[i];
    InitializeValues(i);
}

//__declspec(noinline)
void cFastGen::MarkWindow(cSieveWindow& sieve_window, cSeedElements& seed_elements)
{
    sieve_window.fill(0xFF);
    tpPrime n = firstPrimes[root_N];
    for (auto& element : seed_elements)
    {
        assert(n == element.prime_value);
        tpPrime position = element.next_pos;
        while (position < BITS_PER_WINDOW)
        {
            //mark composite
            sieve_window[position >> 3] &= bit_inverse_masks[position & (8-1)];
                //unsigned byte = position / 8;
                //unsigned bit = position % 8;
                //sieve_window[byte] &= bit_inverse_masks[bit];
            //advance
            position += n;
        }
        element.next_pos = position - BITS_PER_WINDOW;
        n += element.prime_step;
    }
}

//__forceinline
tpPrime cFastGen::CountWindows(cAllWindows& windows, const tpPrime window_start)
{
    tpPrime numPrimes = 0;

    tpPrime n = window_start;
    tpPrime k = 0;
    bool stop_work = false;

    auto ProcessBit = [&](const cSieveWindow& sieve_window, const uint8_t r, const uint8_t bit)
    {
        if (sieve_window[k] & bit_masks[bit])
        {
            tpPrime nk = n + r;
            if (nk < limit)
            {
                assert((nk == 1) or checker.check_next_prime(nk));
                numPrimes++;
            }
            else
                stop_work = true;
        }
    };

    for (; k < BYTES_PER_WINDOW; k++)
    {
        for (uint8_t b = 0; b < 8; b++)
        {
            ProcessBit(windows.window_01, 1, b);
            ProcessBit(windows.window_07, 7, b);
            ProcessBit(windows.window_11, 11, b);
            ProcessBit(windows.window_13, 13, b);
            ProcessBit(windows.window_17, 17, b);
            ProcessBit(windows.window_19, 19, b);
            ProcessBit(windows.window_23, 23, b);
            ProcessBit(windows.window_29, 29, b);

            n += root_primorial;
        }
        if (stop_work)
            break;
    };

    return numPrimes;
}

tpPrime cFastGen::CountWindow(cSieveWindow& sieve_window, const uint8_t r, const tpPrime window_start)
{
    tpPrime numPrimes = 0;

    tpPrime n = window_start;
    tpPrime k = 0;
    bool stop_work = false;

    auto ProcessBit = [&](const uint8_t bit)
    {
        if (sieve_window[k] & bit_masks[bit])
        {
            tpPrime nk = n + r;
            if (nk < limit)
                numPrimes++;
            else
                stop_work = true;
        }
    };

    for (; k < BYTES_PER_WINDOW; k++)
    {
        for (uint8_t b = 0; b < 8; b++)
        {
            ProcessBit(b);

            n += 30;
        }
        if (stop_work)
            break;
    };

    return numPrimes;
}

tpPrime cFastGen::CountWindowsSimple(cAllWindows& windows, const tpPrime window_start, const unsigned win_sz)
{
    tpPrime numPrimes = 0, k = 0;
    for (; k < win_sz; k++)
    {
        numPrimes += 
            bit1_count_table[windows.window_01[k]] +
            bit1_count_table[windows.window_07[k]] +
            bit1_count_table[windows.window_11[k]] +
            bit1_count_table[windows.window_13[k]] +
            bit1_count_table[windows.window_17[k]] +
            bit1_count_table[windows.window_19[k]] +
            bit1_count_table[windows.window_23[k]] +
            bit1_count_table[windows.window_29[k]]            ;
    };

    tpPrime n;
    auto ProcessBit = [&](cSieveWindow& sieve_window, const unsigned char msk, const unsigned char r)
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
    auto ProcessBits = [&](const unsigned char msk)
    {
        ProcessBit(windows.window_01, msk, 1);
        ProcessBit(windows.window_07, msk, 7);
        ProcessBit(windows.window_11, msk, 11);
        ProcessBit(windows.window_13, msk, 13);
        ProcessBit(windows.window_17, msk, 17);
        ProcessBit(windows.window_19, msk, 19);
        ProcessBit(windows.window_23, msk, 23);
        ProcessBit(windows.window_29, msk, 29);
    };
    if (win_sz < BYTES_PER_WINDOW)
    {
        n = window_start + 240 * k;
        ProcessBits(1); n += 30;
        ProcessBits(2); n += 30;
        ProcessBits(4); n += 30;
        ProcessBits(8); n += 30;
        ProcessBits(16); n += 30;
        ProcessBits(32); n += 30;
        ProcessBits(64); n += 30;
        ProcessBits(128);
        //ProcessBit(windows.window_01, 1, 1);
        //ProcessBit(windows.window_07, 1, 7);
        //ProcessBit(windows.window_11, 1, 11);
        //ProcessBit(windows.window_13, 1, 13);
        //ProcessBit(windows.window_17, 1, 17);
        //ProcessBit(windows.window_19, 1, 19);
        //ProcessBit(windows.window_23, 1, 23);
        //ProcessBit(windows.window_29, 1, 29);
        //n += 30;
        //ProcessBit(windows.window_01, 2, 1);
        //ProcessBit(windows.window_07, 2, 7);
        //ProcessBit(windows.window_11, 2, 11);
        //ProcessBit(windows.window_13, 2, 13);
        //ProcessBit(windows.window_17, 2, 17);
        //ProcessBit(windows.window_19, 2, 19);
        //ProcessBit(windows.window_23, 2, 23);
        //ProcessBit(windows.window_29, 2, 29);
        //n += 30;
        //ProcessBit(windows.window_01, 4, 1);
        //ProcessBit(windows.window_07, 4, 7);
        //ProcessBit(windows.window_11, 4, 11);
        //ProcessBit(windows.window_13, 4, 13);
        //ProcessBit(windows.window_17, 4, 17);
        //ProcessBit(windows.window_19, 4, 19);
        //ProcessBit(windows.window_23, 4, 23);
        //ProcessBit(windows.window_29, 4, 29);
        //n += 30;
        //ProcessBit(windows.window_01, 8, 1);
        //ProcessBit(windows.window_07, 8, 7);
        //ProcessBit(windows.window_11, 8, 11);
        //ProcessBit(windows.window_13, 8, 13);
        //ProcessBit(windows.window_17, 8, 17);
        //ProcessBit(windows.window_19, 8, 19);
        //ProcessBit(windows.window_23, 8, 23);
        //ProcessBit(windows.window_29, 8, 29);
        //n += 30;
        //ProcessBit(windows.window_01, 16, 1);
        //ProcessBit(windows.window_07, 16, 7);
        //ProcessBit(windows.window_11, 16, 11);
        //ProcessBit(windows.window_13, 16, 13);
        //ProcessBit(windows.window_17, 16, 17);
        //ProcessBit(windows.window_19, 16, 19);
        //ProcessBit(windows.window_23, 16, 23);
        //ProcessBit(windows.window_29, 16, 29);
        //n += 30;
        //ProcessBit(windows.window_01, 32, 1);
        //ProcessBit(windows.window_07, 32, 7);
        //ProcessBit(windows.window_11, 32, 11);
        //ProcessBit(windows.window_13, 32, 13);
        //ProcessBit(windows.window_17, 32, 17);
        //ProcessBit(windows.window_19, 32, 19);
        //ProcessBit(windows.window_23, 32, 23);
        //ProcessBit(windows.window_29, 32, 29);
        //n += 30;
        //ProcessBit(windows.window_01, 64, 1);
        //ProcessBit(windows.window_07, 64, 7);
        //ProcessBit(windows.window_11, 64, 11);
        //ProcessBit(windows.window_13, 64, 13);
        //ProcessBit(windows.window_17, 64, 17);
        //ProcessBit(windows.window_19, 64, 19);
        //ProcessBit(windows.window_23, 64, 23);
        //ProcessBit(windows.window_29, 64, 29);
        //n += 30;
        //ProcessBit(windows.window_01, 128, 1);
        //ProcessBit(windows.window_07, 128, 7);
        //ProcessBit(windows.window_11, 128, 11);
        //ProcessBit(windows.window_13, 128, 13);
        //ProcessBit(windows.window_17, 128, 17);
        //ProcessBit(windows.window_19, 128, 19);
        //ProcessBit(windows.window_23, 128, 23);
        //ProcessBit(windows.window_29, 128, 29);
    }

    //tpPrime n = window_start;
    //bool stop_work = false;
    //auto ProcessBit = [&](const cSieveWindow& sieve_window, const uint8_t r, const uint8_t bit)
    //{
    //    if (sieve_window[k] & bit_masks[bit])
    //    {
    //        tpPrime nk = n + r;
    //        if (nk < limit)
    //        {
    //            assert((nk == 1) or checker.check_next_prime(nk));
    //            numPrimes++;
    //        }
    //        else
    //            stop_work = true;
    //    }
    //};
    //for (; k < BYTES_PER_WINDOW; k++)
    //{
    //    for (uint8_t b = 0; b < 8; b++)
    //    {
    //        ProcessBit(windows.window_01, 1, b);
    //        ProcessBit(windows.window_07, 7, b);
    //        ProcessBit(windows.window_11, 11, b);
    //        ProcessBit(windows.window_13, 13, b);
    //        ProcessBit(windows.window_17, 17, b);
    //        ProcessBit(windows.window_19, 19, b);
    //        ProcessBit(windows.window_23, 23, b);
    //        ProcessBit(windows.window_29, 29, b);
    //        n += root_primorial;
    //    }
    //    if (stop_work)
    //        break;
    //};

    return numPrimes;
}

tpPrime cFastGen::GenPrimes(const tpPrime lmt, bool dogen)
{
    limit = lmt;

    // 2, 3, 5 are not generated, but 1 is
    tpPrime numPrimes = 3 - 1;

    tpPrime window_start = 0;
    for (; window_start < limit; window_start += WINDOW_SIZE)
    {
        MarkWindow(all_sieve_windows1.window_01, all_seed_elements.seeds_01);
        MarkWindow(all_sieve_windows1.window_07, all_seed_elements.seeds_07);
        MarkWindow(all_sieve_windows1.window_11, all_seed_elements.seeds_11);
        MarkWindow(all_sieve_windows1.window_13, all_seed_elements.seeds_13);
        MarkWindow(all_sieve_windows1.window_17, all_seed_elements.seeds_17);
        MarkWindow(all_sieve_windows1.window_19, all_seed_elements.seeds_19);
        MarkWindow(all_sieve_windows1.window_23, all_seed_elements.seeds_23);
        MarkWindow(all_sieve_windows1.window_29, all_seed_elements.seeds_29);
        numPrimes += dogen ?
            CountWindows(all_sieve_windows1, window_start) :
            CountWindowsSimple(all_sieve_windows1, window_start,
                std::min(BYTES_PER_WINDOW, unsigned((limit - window_start) / 30 / 8)));
    }

    return numPrimes;
}

//__declspec(noinline)
tpPrime cFastGen::GenAllInPrl(const tpPrime lmt)
{
    limit = lmt;

    // 2, 3, 5 are not generated, but 1 is
    num_primes = 3 - 1;

    std::thread tc(&cFastGen::CountTT, this);
    MarkTT();
    tc.join();

    return num_primes;
}

void cFastGen::MarkTT(void)
{
    bool tictoc = true;
    tpPrime window_start_1 = 0;
    tpPrime window_start_2 = WINDOW_SIZE;
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

            MarkWindow(all_sieve_windows1.window_01, all_seed_elements.seeds_01);
            MarkWindow(all_sieve_windows1.window_07, all_seed_elements.seeds_07);
            MarkWindow(all_sieve_windows1.window_11, all_seed_elements.seeds_11);
            MarkWindow(all_sieve_windows1.window_13, all_seed_elements.seeds_13);
            MarkWindow(all_sieve_windows1.window_17, all_seed_elements.seeds_17);
            MarkWindow(all_sieve_windows1.window_19, all_seed_elements.seeds_19);
            MarkWindow(all_sieve_windows1.window_23, all_seed_elements.seeds_23);
            MarkWindow(all_sieve_windows1.window_29, all_seed_elements.seeds_29);
            window_start_1 += 2 * WINDOW_SIZE;

            sync_1.Release();
        }
        else
        {
            if (window_start_2 >= limit)
                break;

            tictoc = true;

            sync_2.WaitTurnAndAcquire(expected_winstate_2);
            expected_winstate_2 += 2;

            MarkWindow(all_sieve_windows2.window_01, all_seed_elements.seeds_01);
            MarkWindow(all_sieve_windows2.window_07, all_seed_elements.seeds_07);
            MarkWindow(all_sieve_windows2.window_11, all_seed_elements.seeds_11);
            MarkWindow(all_sieve_windows2.window_13, all_seed_elements.seeds_13);
            MarkWindow(all_sieve_windows2.window_17, all_seed_elements.seeds_17);
            MarkWindow(all_sieve_windows2.window_19, all_seed_elements.seeds_19);
            MarkWindow(all_sieve_windows2.window_23, all_seed_elements.seeds_23);
            MarkWindow(all_sieve_windows2.window_29, all_seed_elements.seeds_29);
            window_start_2 += 2 * WINDOW_SIZE;

            sync_2.Release();
        }
    } while (true);
}

void cFastGen::CountTT(void)
{
    bool tictoc = true;
    tpPrime window_start_1 = 0;
    tpPrime window_start_2 = WINDOW_SIZE;
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

            num_primes += CountWindowsIterate(all_sieve_windows1, window_start_1);
            window_start_1 += 2 * WINDOW_SIZE;

            sync_1.Release();
        }
        else
        {
            if (window_start_2 >= limit)
                break;

            tictoc = true;

            sync_2.WaitTurnAndAcquire(expected_winstate_2);
            expected_winstate_2 += 2;

            num_primes += CountWindowsIterate(all_sieve_windows2, window_start_2);
            window_start_2 += 2 * WINDOW_SIZE;

            sync_2.Release();
        }
    } while (true);
}

tpPrime cFastGen::CountWindowsIterate(cAllWindows& windows, const tpPrime window_start)
{
    tpPrime numPrimes = 0;

    tpPrime n = window_start;
    tpPrime k = 0;
    bool stop_work = false;

    std::array<cSieveWindow*, 8> wndws = { &windows.window_01, &windows.window_07,
                    &windows.window_11, &windows.window_13, &windows.window_17,
                    &windows.window_19, &windows.window_23,&windows.window_29 };

    for (; k < BYTES_PER_WINDOW; k++)
    {
        for (uint8_t b = 0; b < 8; b++)
        {
            for (int i = 0; i < 8; i++)
                if ((*wndws[i])[k] & bit_masks[b])
                {
                    tpPrime nk = n + root_sieve_values[i];
                    if (nk < limit)
                    {
                        assert((nk == 1) or checker.check_next_prime(nk));
                        numPrimes++;
                    }
                    else
                    {
                        stop_work = true;
                        //break;
                    }
                }

            n += root_primorial;
        }
        if (stop_work)
            break;
    };

    return numPrimes;
}

void cFastGen::MarkNCount(void)
{
    for(tpPrime window_start = 0;;)
    {
        if (window_start >= limit)
            break;
        MarkWindow(all_sieve_windows1.window_01, all_seed_elements.seeds_01);
        MarkWindow(all_sieve_windows1.window_07, all_seed_elements.seeds_07);
        MarkWindow(all_sieve_windows1.window_11, all_seed_elements.seeds_11);
        MarkWindow(all_sieve_windows1.window_13, all_seed_elements.seeds_13);
        MarkWindow(all_sieve_windows1.window_17, all_seed_elements.seeds_17);
        MarkWindow(all_sieve_windows1.window_19, all_seed_elements.seeds_19);
        MarkWindow(all_sieve_windows1.window_23, all_seed_elements.seeds_23);
        MarkWindow(all_sieve_windows1.window_29, all_seed_elements.seeds_29);

        num_primes += CountWindowsIterate(all_sieve_windows1, window_start);

        window_start += WINDOW_SIZE;
    }
}

//__declspec(noinline)
//tpPrime cFastGen::GenAllNotPrl(const tpPrime lmt)
//{
//    limit = lmt;
//
//    // 2, 3, 5 are not generated, but 1 is
//    num_primes = 3 - 1;
//
//    std::thread tc(&cFastGen::CountTT, this);
//    MarkTT();
//    tc.join();
//
//    return num_primes;
//}

void testAdvancedCachedSieve(const uint64_t LIMIT, bool dogen)
{
    if (sqrt(LIMIT) > prms_val_lmt[LMT_USED]) __debugbreak();

    std::cout << "\n - Cached | Advanced "
        << (dogen ? "w.gen.primes" : "count only") << " - \n";
    
    cTimer tmr;
    std::vector<decltype(tmr.LapTime())> times;

    tmr.Start();

    for (auto i : Range<1>())
    {
        cFastGen worker;
        tpPrime n = worker.GenPrimes(LIMIT, dogen);  nln();
        std::cout << n << " primes up to " << LIMIT; nln();
        times.push_back(tmr.LapTime(true));
    }

    tmr.Stop();
    nln(true);
    std::cout << "Average compute time: " << Average(times);
    nln(true);
}

void testAdvancedCachedSieve2T(const uint64_t LIMIT)
{
    if (sqrt(LIMIT) > prms_val_lmt[LMT_USED]) __debugbreak();

    std::cout << "\n - Cached | Advanced 2 threads - \n";

    cTimer tmr;
    std::vector<decltype(tmr.LapTime())> times;

    tmr.Start();

    for (auto i : Range<1>())
    {
        cFastGen worker;
        tpPrime n = worker.GenAllInPrl(LIMIT);  nln();
        std::cout << n << " primes up to " << LIMIT; nln();
        times.push_back(tmr.LapTime(true));
    }

    tmr.Stop();
    nln(true);
    std::cout << "Average compute time: " << Average(times);
    nln(true);
}