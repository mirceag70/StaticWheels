#include "Helper.h"

					//			-- Primorials --
					// P( 2) =							6	 5
					// P( 3) =						   30	 7
					// P( 4) =						  210	11
					// P( 5) =						2,310	13
					// P( 6) =					   30,030	17
					// P( 7) =			          510,510	19
					// P( 8) =					9,699,690	23
					// P( 9) =				  223,092,870	29
					// P(10) =			    6,469,693,230	31
					// P(11) =		      200,560,490,130	37
					// P(12) =		    7,420,738,134,810
					// P(13) =        304,250,263,527,210
					// P(14) =	   13,082,761,331,670,030
					// P(15) =	  614,889,782,588,491,410
					// P(16) = 14,142,414,403,480,493,114

tpPrime primorial;
unsigned char root_sieve[values4N[10] + 1];
tpPrime root_sieve_size = 0;
unsigned root_sieve_N = 0;

//brute force
void InitializeRootSieveFix(tpPrime root_sieve_tmp[], const std::function<bool(tpPrime)>& test)
{
    root_sieve_size = 0;
    unsigned char step = 4;
    for (tpPrime i = 1; i < primorial; i += step, step = 6 - step) 
        if (test(i)) 
            root_sieve_tmp[root_sieve_size++] = i;
}
void InitiliazeRoots(const unsigned char N, tpPrime root_sieve_tmp[])
{
    std::cout << " - brute force";
    switch (N)
    {
    case 3:
        InitializeRootSieveFix(root_sieve_tmp,
            [](tpPrime n) { return (n % 5); });
        break;
    case 4:
        InitializeRootSieveFix(root_sieve_tmp,
            [](tpPrime n) { return (n % 5) and (n % 7); });
        break;
    case 5:
        InitializeRootSieveFix(root_sieve_tmp,
            [](tpPrime n) { return (n % 5) and (n % 7) and (n % 11); });
        break;
    case 6:
        InitializeRootSieveFix(root_sieve_tmp,
            [](tpPrime n) { return (n % 5) and (n % 7) and (n % 11)
            and (n % 13); });
        break;
    case 7:
        InitializeRootSieveFix(root_sieve_tmp,
            [](tpPrime n) { return (n % 5) and (n % 7) and (n % 11)
            and (n % 13) and (n % 17); });
        break;
    case 8:
        InitializeRootSieveFix(root_sieve_tmp,
            [](tpPrime n) { return (n % 5) and (n % 7) and (n % 11)
            and (n % 13) and (n % 17) and (n % 19); });            
        break;
    case 9:
        InitializeRootSieveFix(root_sieve_tmp,
            [](tpPrime n) { return (n % 5) and (n % 7) and (n % 11)
            and (n % 13) and (n % 17) and (n % 19) and (n % 23); });            
            //[](tpPrime n) { return (n>5*(n/5)) and (n > 7 * (n / 7)) and (n > 11 * (n / 11))
            //and (n > 13 * (n / 13)) and (n > 17 * (n / 17)) and (n > 19 * (n / 19)) and (n > 23 * (n / 23)); });
        break;
    case 10:
        InitializeRootSieveFix(root_sieve_tmp,
            [](tpPrime n) { return (n % 5) and (n % 7) and (n % 11)
            and (n % 13) and (n % 17) and (n % 19) and (n % 23) and (n % 29); });
        break;
    default:
        assert(false);
    }
}
void InitializeRootSieve1(const unsigned char N)
{
    if (not((N > 2) and (N <= 10))) return;

    tpPrime i;

    std::cout << "\nInitialize Root Sieve for N = " << (unsigned)N;


    cTimer tmr;
    tmr.Start();

    primorial = 1;
    for (i = 0; i < N; i++) primorial *= firstPrimes[i];

    tpPrime* root_sieve_tmp = new tpPrime[values4N[N]];

    InitiliazeRoots(N, root_sieve_tmp);
    assert(root_sieve_size == values4N[N]);

#pragma warning(push)
#pragma warning(disable:6385)
    // transform absolute values in step increments
    unsigned gap, maxV = 0;
    for (i = 1; i < root_sieve_size; i++)
    {
        gap = (unsigned)(root_sieve_tmp[i] - root_sieve_tmp[i - 1]);
        assert(gap < 256); root_sieve[i - 1] = gap;
        if (root_sieve[i - 1] > maxV) maxV = root_sieve[i - 1];
    }
#pragma warning(disable:6001)
    gap = (unsigned)(primorial + 1 - root_sieve_tmp[root_sieve_size - 1]);
    assert(gap < 256); root_sieve[root_sieve_size - 1] = gap;
    if (root_sieve[root_sieve_size - 1] > maxV) maxV = root_sieve[root_sieve_size - 1];
#pragma warning(pop)

    delete[] root_sieve_tmp;

    tmr.Stop(true, "wheel init");

    std::cout << "\nValues: " << root_sieve_size << " / Primorial: " << primorial;
    std::cout << " / Max gap: " << maxV; nln();
    for (auto i : Range<8>()) std::cout << (unsigned)(root_sieve[i]) << " | ";
}

//6k sieving
void SieveRoots(const unsigned char N, tpPrime root_sieve_tmp[])
{
    std::cout << " - 6k 1bit";

    const tpPrime svlen = primorial / 24 + 1;
    unsigned char* sv = new unsigned char[svlen];
    memset(sv, 0, svlen);

    for (unsigned char i = 2; i < N; i++)
    {
        unsigned char step = 4;
        tpPrime p = firstPrimes[i];
        for (tpPrime n = p; n < primorial; n += p * step , step = 6 - step)
            SetBit(no2idx(n), sv);
    }

    root_sieve_size = 0; root_sieve_tmp[root_sieve_size++] = 1;
    tpPrime n = firstPrimes[N];
    unsigned char step = ((n - 6*(n/6)) == 1) ? 4 : 2;
    for (; n < primorial; n += step, step = 6 - step)
        if (!GetBit(no2idx(n), sv))
            root_sieve_tmp[root_sieve_size++] = n;
            //root_sieve_size++;

    //std::cout << root_sieve_size; nln();

    delete[] sv;
}
void InitializeRootSieve2(const unsigned char N, bool vb = false)
{
    if (not((N > 2) and (N <= 10))) return;

    tpPrime i;
    root_sieve_N = N;

    if (vb) std::cout << "\nInitialize Root Sieve for N = " << (unsigned)N;


    cTimer tmr;
    tmr.Start();

    primorial = 1;
    for (i = 0; i < N; i++) primorial *= firstPrimes[i];

    tpPrime* root_sieve_tmp = new tpPrime[values4N[N]];

    SieveRoots(N, root_sieve_tmp); 
    assert(root_sieve_size == values4N[N]);
    
#pragma warning(push)
#pragma warning(disable:6385)
    // transform absolute values in step increments
    unsigned gap, maxV = 0;
    for (i = 1; i < root_sieve_size; i++)
    {
        gap = (unsigned)(root_sieve_tmp[i] - root_sieve_tmp[i - 1]);
        assert(gap < 256); root_sieve[i - 1] = gap;
        if (root_sieve[i - 1] > maxV) maxV = root_sieve[i - 1];
    }
//#pragma warning(disable:6001)
    gap = (unsigned)(primorial + 1 - root_sieve_tmp[root_sieve_size - 1]);
    assert(gap < 256); root_sieve[root_sieve_size - 1] = gap;
    if (root_sieve[root_sieve_size - 1] > maxV) maxV = root_sieve[root_sieve_size - 1];
#pragma warning(pop)

    delete[] root_sieve_tmp;

    tmr.Stop(true, "wheel init");

    if (vb)
    {
        std::cout << "\nValues: " << root_sieve_size << " / Primorial: " << primorial;
        std::cout << " / Max gap: " << maxV; nln();
        for (auto i : Range<8>()) std::cout << (unsigned)(root_sieve[i]) << " | ";
    }
}

//pritchard
void SoPRoots(const unsigned char N, tpPrime root_sieve_tmp[])
{
    const uint64_t i_max = primorial / 24;
    uint8_t* flags = new uint8_t[i_max + 2];
    uint8_t* s = new uint8_t[primorial];

    std::cout << " - SoP";

    auto SetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        flags[idx] |= BIT_MASK[bit];
    };
    auto ResetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        flags[idx] &= BIT_RESET_MASK[bit];
    };
    auto GetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        return (flags[idx] & BIT_MASK[bit]);
    };

    uint64_t i, p, i_maxS = 0, length = 6;
    for (i = 0; i <= i_max; i++) flags[i] = 0;

    auto Delete = [&](uint64_t pf) { ResetBit(no2idx(pf)); };
    auto next = [&](uint64_t w)
    {
        uint64_t i_w = no2idx(w);
        uint64_t next_i_w = i_w + s[i_w];
        if (!GetBit(next_i_w))
        {
            next_i_w += s[next_i_w];
            s[i_w] = (uint8_t)(next_i_w - i_w);
        }
        return (idx2no(next_i_w));
    };
    auto next_simple = [&](uint64_t w)
    {
        uint64_t i_w = no2idx(w);
        uint64_t next_i_w = i_w + s[i_w];
        return (idx2no(next_i_w));
    };
    auto DeleteMultiples = [&](uint64_t p)
    {
        uint64_t i_f = no2idx(p);
        for (uint64_t f = p; p * f <= length; f = idx2no(i_f))
        {
            Delete(p * f);
            i_f += s[i_f];
        }
    };
    auto Append = [&](uint64_t w)
    {
        uint64_t i_w = no2idx(w);
        SetBit(i_w);    // prime candidate
        s[i_maxS] = (uint8_t)(i_w - i_maxS); i_maxS = i_w;
    };
    //auto Append_simple = [&](uint64_t w)
    //{
    //    uint64_t i_w = no2idx(w);
    //    SetBit(i_w);    // prime candidate
    //};
    auto ExtendToNextP = [&](uint64_t p)
    {
        uint64_t base = length;
        Append(base + 1);
        for (uint64_t w = p; w < length; w = next(w))
            Append(base + w);
        base += length;
        for (int i = 2; i < p; i++)
        {
            Append(base + 1);
            for (uint64_t w = p; w < length; w = next_simple(w))
                Append(base + w);
            base += length;
        }
        length = base;
    };
    //auto ExtendToN = [&](uint64_t p)
    //{
    //    uint64_t base = length;
    //    Append(base + 1);
    //    for (uint64_t w = p; w < length; w = next(w))
    //        if (base + w <= N)
    //            Append_simple(base + w);
    //        else
    //            break;
    //    base += length;
    //    for (int i = 2; i < p and base < N; i++)
    //    {
    //        Append_simple(base + 1);
    //        for (uint64_t w = p; w < length; w = next_simple(w))
    //            if (base + w <= N)
    //                Append_simple(base + w);
    //            else
    //                break;
    //        base += length;
    //    }
    //
    //    length = N;
    //};

    Append(5); 

    for (i = 2; i < N; i++)
    {
        p = firstPrimes[i];
        ExtendToNextP(p);
        DeleteMultiples(p);
    }

    root_sieve_size = 0; root_sieve_tmp[root_sieve_size++] = 1;

    //i = (p - 5) / 24;
    //uint64_t pbase = (i * 24) + 5;
    //uint8_t step{}, b = (uint8_t)(p - pbase);
    //switch (b)
    //{
    //case 0: b = 0; step = 2; break; case 12: b = 4; step = 2; break;
    //case 2: b = 1; step = 4; break; case 14: b = 5; step = 4; break;
    //case 6: b = 2; step = 2; break; case 18: b = 6; step = 2; break;
    //case 8: b = 3; step = 4; break; case 20: b = 7; step = 4; break;
    //default: assert(false);
    //}
    //uint8_t flgs = flags[i];
    //assert(flgs & BIT_MASK[b]);
    //for (; b < 8; b++)
    //{
    //    if (flgs & BIT_MASK[b]) { numPrimes++; AddPrime(p); }
    //    p += step; step = 6 - step;
    //}
    //if (N < 4) { root_sieve_tmp[root_sieve_size++] = 5; }; 
    if (N < 4) { root_sieve_tmp[root_sieve_size++] = 7; };
    if (N < 5) { root_sieve_tmp[root_sieve_size++] = 11; };
    if (N < 6) { root_sieve_tmp[root_sieve_size++] = 13; };
    if (N < 7) { root_sieve_tmp[root_sieve_size++] = 17; };
    if (N < 8) { root_sieve_tmp[root_sieve_size++] = 19; };
    if (N < 9) { root_sieve_tmp[root_sieve_size++] = 23; };
    if (N < 10) { root_sieve_tmp[root_sieve_size++] = 29; };
    uint8_t flgs = flags[1]; p = 29;
    if (flgs & BIT_MASK[1]) { root_sieve_tmp[root_sieve_size++] = p + 2; }
    if (flgs & BIT_MASK[2]) { root_sieve_tmp[root_sieve_size++] = p + 6; }
    if (flgs & BIT_MASK[3]) { root_sieve_tmp[root_sieve_size++] = p + 8; }
    if (flgs & BIT_MASK[4]) { root_sieve_tmp[root_sieve_size++] = p + 12; }
    if (flgs & BIT_MASK[5]) { root_sieve_tmp[root_sieve_size++] = p + 14; }
    if (flgs & BIT_MASK[6]) { root_sieve_tmp[root_sieve_size++] = p + 18; }
    if (flgs & BIT_MASK[7]) { root_sieve_tmp[root_sieve_size++] = p + 20; }
    for (i = 2, p += 24; i < i_max - 1; i++, p += 24)
    {
        uint8_t flgs = flags[i];
        if (flgs & BIT_MASK[0]) { root_sieve_tmp[root_sieve_size++] = p; }
        if (flgs & BIT_MASK[1]) { root_sieve_tmp[root_sieve_size++] = p + 2; }
        if (flgs & BIT_MASK[2]) { root_sieve_tmp[root_sieve_size++] = p + 6; }
        if (flgs & BIT_MASK[3]) { root_sieve_tmp[root_sieve_size++] = p + 8; }
        if (flgs & BIT_MASK[4]) { root_sieve_tmp[root_sieve_size++] = p + 12; }
        if (flgs & BIT_MASK[5]) { root_sieve_tmp[root_sieve_size++] = p + 14; }
        if (flgs & BIT_MASK[6]) { root_sieve_tmp[root_sieve_size++] = p + 18; }
        if (flgs & BIT_MASK[7]) { root_sieve_tmp[root_sieve_size++] = p + 20; }
    }
    while (p <= primorial)
    {
        uint8_t flgs = flags[i++], step = 2;
        for (uint8_t b = 0; b < 8; b++)
        {
            if (p > primorial) break;
            if (flgs & BIT_MASK[b]) { root_sieve_tmp[root_sieve_size++] = p; }
            p += step; step = 6 - step;
        }
    }

    delete[] flags;
    delete[] s;
}
void InitializeRootSieve3(const unsigned char N)
{
    if (not((N > 2) and (N <= 10))) return;

    tpPrime i;

    std::cout << "\nInitialize Root Sieve for N = " << (unsigned)N;


    cTimer tmr;
    tmr.Start();

    primorial = 1;
    for (i = 0; i < N; i++) primorial *= firstPrimes[i];

    tpPrime* root_sieve_tmp = new tpPrime[values4N[N]];

    SoPRoots(N, root_sieve_tmp);
    assert(root_sieve_size == values4N[N]);

#pragma warning(push)
#pragma warning(disable:6385)
    // transform absolute values in step increments
    unsigned gap, maxV = 0;
    for (i = 1; i < root_sieve_size; i++)
    {
        gap = (unsigned)(root_sieve_tmp[i] - root_sieve_tmp[i - 1]);
        assert(gap < 256); root_sieve[i - 1] = gap;
        if (root_sieve[i - 1] > maxV) maxV = root_sieve[i - 1];
    }
#pragma warning(disable:6001)
    gap = (unsigned)(primorial + 1 - root_sieve_tmp[root_sieve_size - 1]);
    assert(gap < 256); root_sieve[root_sieve_size - 1] = gap;
    if (root_sieve[root_sieve_size - 1] > maxV) maxV = root_sieve[root_sieve_size - 1];
#pragma warning(pop)

    delete[] root_sieve_tmp;

    tmr.Stop(true, "wheel init");

    std::cout << "\nValues: " << root_sieve_size << " / Primorial: " << primorial;
    std::cout << " / Max gap: " << maxV; nln();
    for (auto i : Range<8>()) std::cout << (unsigned)(root_sieve[i]) << " | ";
}
