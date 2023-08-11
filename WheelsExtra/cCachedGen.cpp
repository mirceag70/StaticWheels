#include "cCachedGen.h"
//#include "cSeedPrimes.h"

tpPrime cCachedGen::SoE_Rooted_Simple(void)
{
    tpPrime n, n1;
    uint8_t i, j;

    nln(true); std::cout << "8 Rooted SoE";

    sieve = new std::remove_pointer<decltype(sieve)>::type[limit];
    for (tpPrime i = 0; i < limit; i++) sieve[i] = 0xFF;

    unsigned numPrimes = 0;
    //to account for first N primes
    for (i = 0; i < root_N; i++)
    {
        checker.check_next_prime(firstPrimes[i]);
        numPrimes++;
    }

    n = 1ull + root_sieve_steps[0];
    i = 1;
    do
    {
        if (sieve[n])
        {
            assert(checker.check_next_prime(n));
            numPrimes++;

            n1 = n * (1ull + root_sieve_steps[0]);
            j = 1;
            do
            {
                sieve[n1] = false;
                n1 += n * root_sieve_steps[j];
                IncrementRootIdx(j);
            } while (n1 < limit);
        }
        n += root_sieve_steps[i];
        IncrementRootIdx(i);
    } while (n < sqrt(limit));

    cTimer tmr;
    tmr.Start();

    //count the rest of primes
    //continue with previous i and n
    do
    {
        if (sieve[n])
        {
            assert(checker.check_next_prime(n));
            numPrimes++;
        }
        n += root_sieve_steps[i];
        IncrementRootIdx(i);
    } while (n < limit);

    tmr.Stop(true); nln();

    if (sieve)
    {
        delete sieve;
        sieve = nullptr;
    }

    return numPrimes;
}

tpPrime cCachedGen::SoE_Rooted_30(void)
{
    uint8_t i, j;
    tpPrime n, n1;
    sieve_size = limit / 30 + ((limit % 30) > 0);
    nln(true); std::cout << "All In Rooted SoE - 30k";

    sieve = new std::remove_pointer<decltype(sieve)>::type[sieve_size];
    for (tpPrime i = 0; i < sieve_size; i++) sieve[i] = 0xFF;
    //mark 1 and numbers outside limit as composite
    RESETbt(sieve, 1);
    if ((limit % 30) > 0)
        for (i = 0; i < (30 - (limit % 30)); i++)
        {
            switch ((limit + i) % 30)
            {
            case  1: case  7: case 11: case 13:
            case 17: case 19: case 23: case 29:
                RESETbt(sieve, limit + i);
            }
        }

    n = 1ull + root_sieve_steps[0];
    i = 1;
    uint8_t ix = 1;
    do
    {
        if (GETbt(sieve,n))
        {
            //n1 = n * (1ull + root_sieve_steps[0]);
            //j = 1;
            n1 = n * n;
            while (root_sieve_values[ix] != n % 30)
                if (++ix == 8)
                    ix = 0;
            j = ix;
            do
            {
                RESETbt(sieve, n1);
                if (n1 == 2093)
                    n = n;
#pragma warning(suppress:28020) // expression not true at this call
                n1 += n * root_sieve_steps[j];
                IncrementRootIdx(j);
            } while (n1 < limit);
        }
        n += root_sieve_steps[i];
        IncrementRootIdx(i);
    } while (n < sqrt(limit));

    tpPrime numPrimes = 0;
    //to account for first N primes
    for (i = 0; i < root_N; i++)
    {
        checker.check_next_prime(firstPrimes[i]);
        numPrimes++;
    }

    numPrimes += CountSieve2();

    if (sieve)
    {
        delete[] sieve;
        sieve = nullptr;
    }

    return numPrimes;
}

tpPrime cCachedGen::CountSieve2(void)
{
    unsigned numPrimes = 0;
    tpPrime n = 0;
    tpPrime k = 0;
    tpPrime _sum = 0;

    tpPrime last_n = 7;
    //unsigned max_diff = 0;

    auto ProcessBit = [&](const unsigned char msk, const unsigned char r)
    {
        if (sieve[k] & msk)
        {
            //if ((n + r) < limit)
            {
                tpPrime nk = n + r;
                _sum += nk;
                //max_diff = std::max(max_diff, unsigned(nk - last_n)); last_n = nk;
                assert(checker.check_next_prime(nk));
                numPrimes++;
            }
        }
    };

    cTimer tmr;
    tmr.Start();

    for (; k < sieve_size; k++)
    {
        ProcessBit(  1,  1); ProcessBit(  2,  7); 
        ProcessBit(  4, 11); ProcessBit(  8, 13); 
        ProcessBit( 16, 17); ProcessBit( 32, 19);
        ProcessBit( 64, 23); ProcessBit(128, 29);
        
        n += 30;
    };

    tmr.Stop(true); nln();
    //std::cout << max_diff; 

    return numPrimes;
}

tpPrime cCachedGen::CountSieve1(void)
{
    unsigned numPrimes = 0;
    tpPrime n = 0;


    cTimer tmr;
    tmr.Start();

    for (tpPrime k = 0; k < sieve_size; k++)
    {
        numPrimes += UnpackSlot(sieve[k], n, 0, 1, 2, 3, 4, 5, 6, 7);
        n += 30;
    };

    tmr.Stop(true); nln();

    return numPrimes;
}

template<typename... Args>
unsigned cCachedGen::UnpackSlot(const tpSieveElement slot, tpPrime nk, const Args&&... rargs)
{
    auto ProcessBit = [&](const unsigned char r) -> unsigned
    {
        if (slot & v_msk[r])
        {
            tpPrime n = nk + root_sieve_values[r];
            if (n > 2)//< limit)
            {
                assert(checker.check_next_prime(n));
                return 1u;
            }
        }
        return 0u;
    };

    return (ProcessBit(rargs) + ...);
}

template<tpPrime limit>
tpPrime Rooted_30(void)
{
    //add 1 slot if required
    constexpr tpPrime sieve_size = limit / 30 + ((limit % 30) > 0);

    tpPrime n, n1;

    //nln(true); std::cout << "All In Rooted SoE - 30k";

    tpSieveElement* sieve = new tpSieveElement[sieve_size];
    for (tpPrime i = 0; i < sieve_size; i++) sieve[i] = 0xFF;
    //mark 1 and numbers outside limit as composite
    RESETbt(sieve, 1);
    for (tpPrime i = 0; i < limit % 30; i++) RESETbt(sieve, limit + i);

    uint8_t i, j;
    n = 1ull + root_sieve_steps[0];
    i = 1;
    uint8_t ix = 1;
    do
    {
        if (GETbt(sieve, n))
        {
            //n1 = n * (1ull + root_sieve_steps[0]);
            //j = 1;
            n1 = n * n;
            while (root_sieve_values[ix] != n % 30)
                if (++ix == 8)
                    ix = 0;
            j = ix;
            do
            {
                RESETbt(sieve, n1);
                n1 += n * root_sieve_steps[j];
                IncrementRootIdx(j);
            } while (n1 < limit);
        }
        n += root_sieve_steps[i];
        IncrementRootIdx(i);
    } while (n < sqrt(limit));

    cChecker checker;

    tpPrime numPrimes = 0;
    //to account for first N primes
    for (i = 0; i < root_N; i++)
    {
        checker.check_next_prime(firstPrimes[i]);
        numPrimes++;
    }

    n = 0;
    tpPrime k = 0;
    tpPrime _sum = 0;

    auto ProcessBit = [&](const unsigned char msk, const unsigned char r)
    {
        if (sieve[k] & msk)
        {
            //if ((n + r) < limit)
            {
                tpPrime nk = n + r;
                _sum += nk;
                //max_diff = std::max(max_diff, unsigned(nk - last_n)); last_n = nk;
                assert(checker.check_next_prime(nk));
                numPrimes++;
            }
        }
    };

    cTimer tmr;
    tmr.Start();

    for (; k < sieve_size; k++)
    {
        ProcessBit(1, 1); ProcessBit(2, 7);
        ProcessBit(4, 11); ProcessBit(8, 13);
        ProcessBit(16, 17); ProcessBit(32, 19);
        ProcessBit(64, 23); ProcessBit(128, 29);

        n += 30;
    };

    if (sieve)
    {
        delete sieve;
        sieve = nullptr;
    }

    return numPrimes;
}
