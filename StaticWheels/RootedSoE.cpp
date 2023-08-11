#include "Helper.h"

uint64_t SoE_6k(const uint64_t Nmax, uint8_t vPrimes[], void*, void*)
{
    auto idx2no = [](uint64_t idx) { return 3 * idx + 5 - (idx & 1); };
    auto no2idx = [](uint64_t no) { return no / 3 - 1; };
    auto SetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        vPrimes[idx] |= BIT_MASK[bit];
    };
    auto GetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        return (vPrimes[idx] & BIT_MASK[bit]);
    };

    const uint64_t i_max = Nmax / (3ull * 8);
    const uint64_t Nsqrt = (uint64_t)floor(sqrt(Nmax));
    unsigned numPrimes = 2; AddPrime(2); AddPrime(3);

    uint64_t i, j;
    for (i = 0; i <= i_max; i++) vPrimes[i] = 0;

    uint8_t stepi = 2;
    for (i = 5; i <= Nsqrt; i += stepi, stepi = 6 - stepi)
        if (!GetBit(no2idx(i)))
        {
            numPrimes++; AddPrime(i);

            uint64_t stepj2 = 2 * i, stepj6 = 6 * i;
            uint64_t stepj = ((i - (i / 3) * 3) == 1) ? (stepj6 - stepj2) : stepj2;

            for (j = i * i; j <= Nmax; j += stepj, stepj = stepj6 - stepj)
                SetBit(no2idx(j));
        }

    for ( /*continue counting primes*/; i <= Nmax; i += stepi, stepi = 6 - stepi)
        if (!GetBit(no2idx(i)))
        {
            numPrimes++; AddPrime(i);
        }

    return numPrimes;
}

void InitializeRootSieve2(const unsigned char N, bool vb = false);
void InitializeRootSieveP(const unsigned char N, bool vb = false);

extern unsigned root_N;
extern unsigned char root_sieve[];
extern tpPrime root_sieve_size;// , root_sieve_N;

unsigned nextWheelIdx(unsigned idx) { return ((idx < (root_sieve_size - 1)) ? idx + 1 : 0); }

uint64_t SoE_rooted(const uint64_t Nmax, uint8_t vPrimes[], void*, void*)
{
    InitializeRootSieve2(root_N);

    //auto idx2no = [](uint64_t idx) { return 3 * idx + 5 - (idx & 1); };
    //auto no2idx = [](uint64_t no) { return no / 3 - 1; };
    auto SetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        vPrimes[idx] |= BIT_MASK[bit];
    };
    auto GetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        return (vPrimes[idx] & BIT_MASK[bit]);
    };

    const uint64_t i_max = Nmax / (3ull * 8);
    const uint64_t Nsqrt = (uint64_t)floor(sqrt(Nmax));
    unsigned numPrimes = root_N;
    for (unsigned i = 0; i < root_N; i++) AddPrime(firstPrimes[i]);

    for (uint64_t i = 0; i <= i_max; i++) vPrimes[i] = 0;

    uint64_t p = 1ull + root_sieve[0];
    unsigned ip = 1;
    for (; p <= Nsqrt; p += root_sieve[ip], ip = nextWheelIdx(ip))
        if (!GetBit(no2idx(p)))
        {
            numPrimes++; AddPrime(p);
            uint64_t n = 1ull + root_sieve[0];
            unsigned in = 1;
            for (; n*p <= Nmax; n += root_sieve[in], in = nextWheelIdx(in))
                SetBit(no2idx(n*p));
        }

    for (/*continue counting primes*/; p <= Nmax; p += root_sieve[ip], ip = nextWheelIdx(ip))
        if (!GetBit(no2idx(p)))
        {
            numPrimes++; AddPrime(p);
        }

    return numPrimes;
}

uint64_t SoE_rooted1(const uint64_t Nmax, uint8_t vPrimes[], void*, void*)
{
    if(root_N <= 7)
        InitializeRootSieve2(root_N);
    else
        InitializeRootSieveP(root_N);

    auto SetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        vPrimes[idx] |= BIT_MASK[bit];
    };
    auto GetBit = [&](uint64_t bitidx)
    {
        uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
        return (vPrimes[idx] & BIT_MASK[bit]);
    };

    const uint64_t i_max = Nmax / (3ull * 8);
    const uint64_t Nsqrt = (uint64_t)floor(sqrt(Nmax));
    for (uint64_t i = 0; i <= i_max; i++) vPrimes[i] = 0;
    unsigned numPrimes = root_N;
    for (unsigned i = 0; i < root_N; i++) AddPrime(firstPrimes[i]);

    uint64_t p = 1ull + root_sieve[0];
    unsigned ip = 1;
    for (; p <= Nsqrt; p += root_sieve[ip], ip = nextWheelIdx(ip))
        if (!GetBit(no2idx(p)))
        {
            numPrimes++; AddPrime(p);
            uint64_t n = p;
            for (unsigned in = ip; n * p <= Nmax; n += root_sieve[in], in = nextWheelIdx(in))
                SetBit(no2idx(n * p));
        }

    for (/*continue counting primes*/; p <= Nmax; p += root_sieve[ip], ip = nextWheelIdx(ip))
        if (!GetBit(no2idx(p)))
        {
            numPrimes++; AddPrime(p);
        }

    return numPrimes;
}

extern tpPrime primorial;
extern unsigned root_sieve_N;

constexpr auto chunk_bytes = 1u << 16; // 2^15 = 32k, cache friendly
constexpr auto chunk_span = chunk_bytes * 4 * 6;
constexpr auto max_chunks = 4200;
unsigned chunks, vals_per_chunk, root_sieve_tmp_sz[max_chunks];
tpPrime* root_sieve_tmpx;// , * root_sieve_tmp;

//6k parallel sieving
void SieveRootsP(const unsigned char N)
{
    const unsigned minsvlen = unsigned(primorial / 24) + 1;
    chunks = minsvlen / chunk_bytes + 1; assert(chunks < max_chunks);
    const unsigned svlen = chunks * chunk_bytes;
    unsigned char* sv = new unsigned char[svlen];
    memset(sv, 0, svlen);

    vals_per_chunk = (unsigned)(1.2 * values4N[N] / chunks);
    //root_sieve_tmp = new tpPrime[vals_per_chunk * chunks];
    root_sieve_tmpx = new tpPrime[vals_per_chunk * chunks];

    //unsigned char* svx = new unsigned char[svlen];
    //memset(svx, 0, svlen);
    //for (unsigned char i = 2; i < N; i++)
    //{
    //    unsigned char step = 4;
    //    tpPrime p = firstPrimes[i];
    //    for (tpPrime n = p; n < primorial; n += p * step, step = 6 - step)
    //        SetBit(no2idx(n), svx);
    //} 

    std::for_each(std::execution::par_unseq,
        range(chunks).begin(), range(chunks).end(), [&](const unsigned chunk)
        {
            const uint64_t start = 5 + uint64_t(chunk) * chunk_span;
            const uint64_t stop = start + chunk_span;

            for (unsigned char i = 2; i < N; i++)
            {
                tpPrime p = firstPrimes[i];
                tpPrime k = 6 * (start / p / 6) + 1;
                tpPrime n = k * p;
                unsigned char step = (n < start) ? (n += 4 * p, 2) : 4;
                //unsigned char step;
                //if (n < start)
                //{
                //    step = 2; n += 4 * p;
                //}
                //else
                //    step = 4;
                for (; n < stop; n += p * step, step = 6 - step)
                    SetBit(no2idx(n), sv);
            }

            const uint64_t idx0 = uint64_t(chunk) * vals_per_chunk;
            tpPrime n = 6 * (start / 6) + 1;
            unsigned char step = (n < start) ? (n += 4, 2) : 4;
            unsigned ix = 0;
            for (; n < stop; n += step, step = 6 - step)
            {
                assert(idx0 + ix < vals_per_chunk * chunks);
                if (!GetBit(no2idx(n), sv))
                    root_sieve_tmpx[idx0 + ix++] = n;
            }
            root_sieve_tmp_sz[chunk] = ix;
        });

    //root_sieve_size = 0; root_sieve_tmp[root_sieve_size++] = 1;
    //tpPrime n = firstPrimes[N];
    //unsigned char step = ((n - 6 * (n / 6)) == 1) ? 4 : 2;
    //for (; n < primorial; n += step, step = 6 - step)
    //{
    //    //assert(GetBit(no2idx(n), sv) == GetBit(no2idx(n), svx));
    //    if (!GetBit(no2idx(n), sv))
    //        root_sieve_tmp[root_sieve_size++] = n;
    //}
    //std::cout << root_sieve_size; nln();
    //tmr.LapTime(true, "prl wheel counting"); 

    delete[] sv;
}
void InitializeRootSieveP(const unsigned char N, bool vb)
{
    if (not((N >= 8) and (N <= 10))) return;
    root_sieve_N = N;

    if (vb) std::cout << "\nInitialize Root Sieve for N = " << (unsigned)N;
    cTimer tmr; tmr.Start();

    primorial = 1;
    for (tpPrime i = 0; i < N; i++) primorial *= firstPrimes[i];

    SieveRootsP(N); tmr.LapTime(true, "sieving and counting");

    // transform absolute values in step increments
    std::for_each(std::execution::par_unseq,
        range(chunks).begin(), range(chunks).end(), [&](const unsigned chunk)
        {
            unsigned idx = 0;
            for (unsigned k = 0; k < chunk; k++) idx += root_sieve_tmp_sz[k];
            const uint64_t idx0 = uint64_t(chunk) * vals_per_chunk;
            uint64_t lastp;
            if (chunk == 0)
                lastp = 1;
            else
            {
                uint64_t idxlastp = idx0 - vals_per_chunk + (uint64_t)root_sieve_tmp_sz[chunk - 1];
                //std::cout << idxlastp << "\n";
                lastp = root_sieve_tmpx[idxlastp - 1];
            }
            for (unsigned ix = 0; ix < root_sieve_tmp_sz[chunk] and idx <= values4N[N]; ix++)
            {
                unsigned gap = unsigned(root_sieve_tmpx[idx0 + ix] - lastp);
                assert(gap < 256); root_sieve[idx++] = gap;
                lastp = root_sieve_tmpx[idx0 + ix];
            }
        });

    root_sieve_size = values4N[N];
    //unsigned gap, idx = 0;
    //uint64_t lastp = 1;
    //for (unsigned chunk = 0; chunk < chunks; chunk++)
    //{
    //    const uint64_t idx0 = uint64_t(chunk) * vals_per_chunk;
    //    for (unsigned ix = 0; ix < root_sieve_tmp_sz[chunk]; ix++)
    //    {
    //        unsigned gap = root_sieve_tmpx[idx0 + ix] - lastp;
    //        assert(root_sieve[idx++] == gap);
    //        lastp = root_sieve_tmpx[idx0 + ix];
    //    }
    //}

    //check values
//    for (tpPrime i = 1; i < root_sieve_size; i++)
//    {
//#pragma warning(push)
//#pragma warning(disable:6385)
//        gap = (unsigned)(root_sieve_tmp[i] - root_sieve_tmp[i - 1]);
//#pragma warning(pop)
//        assert(gap < 256); //root_sieve[i - 1] = gap;
//        assert(root_sieve[i - 1] == gap);
//        if (root_sieve[i - 1] > maxV) maxV = root_sieve[i - 1];
//    }
    //gap = (unsigned)(primorial + 1 - root_sieve_tmp[root_sieve_size - 1]);
    //assert(gap < 256); root_sieve[root_sieve_size - 1] = gap;
    //if (root_sieve[root_sieve_size - 1] > maxV) maxV = root_sieve[root_sieve_size - 1];

    delete[] root_sieve_tmpx;

    tmr.Stop(true, "6k 1bit root sieve");

    if (vb)
    {
        std::cout << "\nValues: " << root_sieve_size << " / Primorial: " << primorial;
        for (auto i : Range<8>()) std::cout << (unsigned)(root_sieve[i]) << " | ";
    }
}