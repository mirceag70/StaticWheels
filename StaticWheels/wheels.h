#pragma once

#include <syncstream>

namespace wheels
{

    template<class NUMBER>
    constexpr NUMBER ConstAbs(NUMBER a)
    {
        static_assert(std::is_arithmetic<NUMBER>::value, "Numeric type required for ConstAbs.");

        return (a >= 0) ? a : (-a);
    }
    template<class NUMBER>
    constexpr long long ConstRound(NUMBER a)
    {
        static_assert(std::is_arithmetic<NUMBER>::value, "Numeric type required for ConstRound.");

        long long rnd = static_cast<long long>(a);
        double delta = static_cast<double>(ConstAbs(a - rnd));

        if (delta <= 0.50)
            return rnd;
        else
            return rnd > 0 ? rnd + 1 : rnd - 1;
    }
    // square root of an integral number using Newtons method
    template <typename INTEGRAL>
    consteval INTEGRAL ConstSqrt(INTEGRAL n)
    {
        static_assert(std::is_integral<INTEGRAL>::value,
            "Integral type required for ConstSqrt.");

        double x = n / 10.0;			// current value
        double gss = (x + (n / x)) / 2;	// current guess 

        while (ConstAbs(gss - x) > 0.25)
            // we need an integer, so our tolerance is quite high
        {
            x = gss;
            gss = (x + (n / x)) / 2;	// better guess
        }

        return static_cast<INTEGRAL>(ConstRound(gss));
    }

    constexpr unsigned _N = 97; static_assert(_N == 97, "_N must be 97 = PI(2^9 i.e. 512)");
    typedef std::array<unsigned, _N> tpFirstPrimesArray;
    auto consteval GenFirst_NPrimes(void)
    {
        tpFirstPrimesArray arr{};

        constexpr unsigned sieve_size = 32;
        unsigned char sieve[sieve_size]{};
        constexpr unsigned max_i = 2ull * sieve_size * sizeof(sieve[0]) * 8;
        for (auto i : range(sieve_size)) sieve[i] = 0;

        for (auto bit : range(1, ConstSqrt(max_i)))
            if (!GetBit(bit, sieve))
            {
                unsigned i = 2 * bit + 1;
                for (unsigned j = i * i; j < max_i; j += 2 * i)
                    SetBit(j / 2, sieve);
            }

        unsigned idx = 1; arr[0] = 2;
        for (auto i : range(1, max_i / 2))
            if (!GetBit(i, sieve))
                arr[idx++] = 2 * i + 1;

        return arr;
    }
    constexpr tpFirstPrimesArray first97primes = GenFirst_NPrimes();

    consteval unsigned Primorial(const unsigned N)
    {
        unsigned prim = 1;
        for (auto i : range(N)) prim *= first97primes[i];
        return prim;
    }

    template<const unsigned _root_N>
    class RootSieve
    {
        static_assert(_root_N < _N, "Primorial() will fail");
        static_assert(_root_N >= 3 and _root_N <= 6, "Bad wheel N");

        static constexpr uint16_t wheel_entries[] = {   8,  48,  480, 5'760 }; 
                                                  //{ 8*1, 8*6, 8*60, 8*720 }
        static bool factors_initialized;

    public:
        static constexpr bool small_wheel = (_root_N <= 4);
        typedef std::conditional<small_wheel, uint8_t, uint16_t>::type tpWheelEntry;

        static constexpr unsigned root_N = _root_N;
        static constexpr unsigned wheel_span = Primorial(root_N);
        static constexpr unsigned wheel_values = wheel_entries[root_N - 3];


        typedef std::array<tpWheelEntry, wheel_values> tpWheelArray;
        static auto consteval GenWheel(void)
        {
            unsigned sieve[wheel_span]{};
            for (auto i : range(wheel_span)) sieve[i] = 1;

            tpWheelArray arr{};

            for (auto i : range(root_N))
            {
                unsigned p = first97primes[i];
                for (unsigned j = p; j < wheel_span; j += p)
                    sieve[j] = 0;
            }

            unsigned idx = 0;
            for (auto i : range(1, wheel_span))
                if (sieve[i]) arr[idx++] = static_cast<tpWheelEntry>(i);
            return arr;
        }
        static constexpr tpWheelArray wheel = GenWheel();

        typedef std::array<tpWheelArray, wheel_values> tpWheelMatrix;
        static tpWheelMatrix mods, factor;

        void GenMods(void)
        {
            std::for_each(std::execution::par_unseq,
                range(wheel_values).begin(), range(wheel_values).end(), [&](const unsigned i)
                {
                    //      Full
                    for (auto j : range(wheel_values))
                        mods[i][j] = (wheel[i] * wheel[j]) % wheel_span;
                    // !due to cache thrashing and concurrency 
                    // !these variants are actually much worse
                    //      Supradiagonal
                    //for (auto j : range(i, wheel_values))
                    //    mods[i][j] = (wheel[i] * wheel[j]) % wheel_span;
                    //      Subdiagonal
                    //for (auto j : range(i))
                    //    mods[i][j] = (wheel[i] * wheel[j]) % wheel_span;
                //std::thread::id t_id = std::this_thread::get_id();
                //std::osyncstream(std::cout) << t_id << "\n";
                });
            //std::for_each(std::execution::par_unseq,
            //    range(wheel_values).begin(), range(wheel_values).end(), [&](const unsigned i)
            //    {
            //        if (i == 0) return;
            //        //      Supradiagonal
            //        //for (auto j : range(i, wheel_values))
            //        //    mods[i][j] = (wheel[i] * wheel[j]) % wheel_span;
            //        
            //        //      Subdiagonal
            //        for (auto j : range(i-1))
            //            mods[i][j] = mods[j][i];
            //   });            
            //for (auto i : range(wheel_values))
            //    for (auto j : range(wheel_values))
            //        mods[i][j] = (wheel[i] * wheel[j]) % wheel_span;
        }

        void GenFactors(void)
        {
            if (factors_initialized) return;
            cTimer tmr; tmr.Start();
            GenMods();
            tmr.LapTime(true, "mods");
            //unsigned values = 0;

            //auto compute_linex = [&](const unsigned i)
            //{
            //    //std::thread::id t_id = std::this_thread::get_id();
            //    //std::osyncstream(std::cout) << t_id << "|";
            //    for (auto j : range(wheel_values))
            //        for (auto k : range(wheel_values))
            //            if (mods[i][k] == wheel[j])
            //            {
            //                factor[i][j] = wheel[k];
            //                break;
            //            }
            //};
            //std::for_each(std::execution::par_unseq,
            //    range(wheel_values).begin(), range(wheel_values).end(), compute_linex);

            auto compute_line = [&](const unsigned i)
            {
                for (auto j : range(wheel_values))
                    if(mods[i][j] == 1)
                    {
                        for (auto k : range(wheel_values))
if constexpr (_root_N == 6)
                            //faster to compute per line and transpose 
                            factor[i][k] = mods[j][k];
else
                            factor[k][i] = mods[j][k];
                        break;
                    }
            };
            std::for_each(std::execution::par_unseq,
                range(wheel_values).begin(), range(wheel_values).end(), compute_line);

            if constexpr (_root_N == 6)
            {
                // transpose for better access pattern in InitLastPK
                constexpr unsigned block_size = (root_N == 6) ? 64 : 8;
                assert(wheel_values % block_size == 0);
                std::for_each(std::execution::par_unseq,
                    range(wheel_values / block_size).begin(), range(wheel_values / block_size).end(), [&](const unsigned iib)
                    {   //cache friendly transposition
                        const unsigned ii = iib * block_size;
                        //transpose the diagonal block
                        for (unsigned i = ii; i < ii + block_size; i++)
                            for (unsigned j = i + 1; j < ii + block_size; j++)
                                std::swap(factor[i][j], factor[j][i]);
                        //transpose the rest of blocks
                        for (unsigned jj = ii + block_size; jj < wheel_values; jj += block_size)
                            for (unsigned i = ii; i < ii + block_size; i++)
                                for (unsigned j = jj; j < jj + block_size; j++)
                                    std::swap(factor[i][j], factor[j][i]);
                    });
            }
            //for (auto i : range(wheel_values))
            //{
                //fprintf(stdout, "%4.1f%%\b\b\b\b\b", i * 100.0 / wheel_values);
                //compute_line(i);
                //unsigned ni = wheel[i];
                //for (auto j : range(wheel_values))
                //{
                //    unsigned nj = wheel[j];
                //    for (auto k : range(wheel_values))
                //    {
                //        if (mods[i][k] == nj)
                //        {
                //            factor[i][j] = wheel[k];
                //            values++; break;
                //        }
                //    }
                //}
            //}
            //assert(values == wheel_values * wheel_values);
            factors_initialized = true;
            tmr.LapTime(true, "factors"); tmr.Stop(true, "full init");
        }

        RootSieve() { GenFactors(); };
    };

    class BigWheel7
    {
        static constexpr uint32_t wheel_entries[] = { 92'160, 1'658'880, 36'495'360, 1'058'365'440 };
        static bool factors_initialized;

    public:
        static constexpr unsigned root_N = 7; static_assert(root_N == 7, "Bad wheel N");
        static constexpr unsigned wheel_span = Primorial(root_N);
        static constexpr unsigned wheel_values = wheel_entries[root_N - 7];

        typedef uint32_t tpWheelEntry;
        typedef std::array<tpWheelEntry, wheel_values> tpWheelArray;
        static tpWheelArray wheel;

        void GenWheel(void)
        {
            const unsigned svlen = wheel_span / 24 + 1;
            unsigned char* sv = new unsigned char[svlen];
            memset(sv, 0, svlen);
            for (unsigned char i = 2; i < root_N; i++)
            {
                unsigned char step = 4;
                tpPrime p = first97primes[i];
                for (tpPrime n = p; n < wheel_span; n += p * step, step = 6 - step)
                    SetBit(no2idx(n), sv);
            }
            unsigned idx0 = 0; wheel[idx0++] = 1;
            tpWheelEntry n = first97primes[root_N];
            unsigned char step = ((n - 6 * (n / 6)) == 1) ? 4 : 2;
            for (; n < wheel_span; n += step, step = 6 - step)
                if (!GetBit(no2idx(n), sv))
                    wheel[idx0++] = n;
            assert(idx0 == wheel_values);
            delete[] sv;
        };
    
        typedef std::array<tpWheelArray*, wheel_values> tpWheelMatrix;
        static tpWheelMatrix mods, factor;
        void GenMods(void)
        {
            // compute supra-diagonal - not anymore...
            std::for_each(std::execution::par_unseq,
                range(wheel_values).begin(), range(wheel_values).end(), [&](const unsigned i)
                {
                    for (auto j : range(wheel_values))
                    {
                        uint64_t nn = wheel[i]; nn *= wheel[j];
                        (*mods[i])[j] = nn % wheel_span;
                    }
                });

            // cache friendly fill up sub-diagonal
            //constexpr unsigned block_size = 64; assert(wheel_values % block_size == 0);
            //std::for_each(std::execution::par_unseq,
            //    range(wheel_values / block_size).begin(), range(wheel_values / block_size).end(), [&](const unsigned iib)
            //    {   
            //        const unsigned ii = iib * block_size;
            //        // the diagonal block
            //        for (unsigned i = ii; i < ii + block_size; i++)
            //            for (unsigned j = i + 1; j < ii + block_size; j++)
            //                (*mods[i])[j] = (*mods[j])[i];
            //        // the rest of blocks
            //        for (unsigned jj = ii + block_size; jj < wheel_values; jj += block_size)
            //            for (unsigned i = ii; i < ii + block_size; i++)
            //                for (unsigned j = jj; j < jj + block_size; j++)
            //                    (*mods[i])[j] = (*mods[j])[i];
            //    });
        }

        //typedef std::array<tpWheelArray*, wheel_values> tpWheelMatrix;
        //static tpWheelMatrix mods, factor;
        //void GenMods(void)
        //{
        //    std::for_each(std::execution::par_unseq,
        //        range(wheel_values).begin(), range(wheel_values).end(), [&](const unsigned i)
        //        {
        //            mods[i] = new tpWheelArray;
        //            for (auto j : range(wheel_values))
        //            {
        //                uint64_t nn = wheel[i]; nn *= wheel[j];
        //                (*mods[i])[j] = nn % wheel_span;
        //            }
        //        });
        //}

        void GenFactors(void)
        {
            if (factors_initialized) return;
            cTimer tmr; tmr.Start();
            GenMods(); tmr.LapTime(true, "mods7");

            auto compute_line = [&](const unsigned i)
            {
                for (auto j : range(wheel_values))
                    if ((*mods[i])[j] == 1) { factor[i] = mods[j]; break; }
                assert(factor[i] != NULL);
            };
            std::for_each(std::execution::par_unseq,
                range(wheel_values).begin(), range(wheel_values).end(), compute_line);
            tmr.LapTime(true, "factors7");

            // transpose for better access pattern in InitLastPK7
            constexpr unsigned block_size = 64; assert(wheel_values % block_size == 0);
            std::for_each(std::execution::par_unseq,
                range(wheel_values / block_size).begin(), range(wheel_values / block_size).end(), [&](const unsigned iib)
                {   //cache friendly transposition
                    const unsigned ii = iib * block_size;
                    //transpose the diagonal block
                    for (unsigned i = ii; i < ii + block_size; i++)
                        for (unsigned j = i + 1; j < ii + block_size; j++)
                            std::swap((*factor[i])[j], (*factor[j])[i]);
                    //transpose the rest of blocks
                    for (unsigned jj = ii + block_size; jj < wheel_values; jj += block_size)
                        for (unsigned i = ii; i < ii + block_size; i++)
                            for (unsigned j = jj; j < jj + block_size; j++)
                                std::swap((*factor[i])[j], (*factor[j])[i]);
                });
            //std::for_each(std::execution::par_unseq,
            //    range(wheel_values).begin(), range(wheel_values).end(), [&](const unsigned i)
            //    {
            //        for (auto j : range(i + 1, wheel_values))
            //        {
            //            std::swap((*factor[i])[j], (*factor[j])[i]);
            //        }
            //    });
            factors_initialized = true;
            tmr.LapTime(true, "transpose7");
            tmr.Stop(true, "full init7"); nln();
        }
 
        static constexpr auto SZ = wheel_values;
        static constexpr auto BLKSZ = 256u; static_assert(SZ% BLKSZ == 0);

        BigWheel7() 
        { 
            cTimer tmr; tmr.Start();
            GenWheel(); tmr.LapTime(true, "wheel7");

            std::for_each(std::execution::par_unseq,
                range(SZ / BLKSZ).begin(), range(SZ / BLKSZ).end(), [&](const unsigned i)
                {
                    tpWheelArray* pv = new tpWheelArray[BLKSZ];
                    for (unsigned j = 0; j < BLKSZ; j++)
                    {
                        assert(uint64_t(i) * BLKSZ + j < UINT32_MAX);
                        mods[i * BLKSZ + j] = pv + j;
                    }
                });

            tmr.Stop(true, "mem7");
        };
        ~BigWheel7() 
        { 
            //cTimer tmr; tmr.Start();

            std::for_each(//std::execution::par_unseq, //FASTER SINGLE THREADED!
                range(SZ / BLKSZ).begin(), range(SZ / BLKSZ).end(), [&](const unsigned i)
                {
                    delete mods[i * BLKSZ];
                });

            //tmr.Stop(true, "free mem7");
        };
    };

//#ifdef TESTWHEELSINITS
    inline bool RootSieve<3>::factors_initialized = false;
    inline    RootSieve<3>::tpWheelMatrix RootSieve<3>::factor;
    inline RootSieve<3>::tpWheelMatrix RootSieve<3>::mods;
    inline bool RootSieve<4>::factors_initialized = false;
    inline RootSieve<4>::tpWheelMatrix RootSieve<4>::factor;
    inline RootSieve<4>::tpWheelMatrix RootSieve<4>::mods;
    inline bool RootSieve<5>::factors_initialized = false;
    inline RootSieve<5>::tpWheelMatrix RootSieve<5>::factor;
    inline RootSieve<5>::tpWheelMatrix RootSieve<5>::mods;
    inline bool RootSieve<6>::factors_initialized = false;
    inline RootSieve<6>::tpWheelMatrix RootSieve<6>::factor;
    inline RootSieve<6>::tpWheelMatrix RootSieve<6>::mods;
//#endif

    inline bool BigWheel7::factors_initialized = false;
    inline BigWheel7::tpWheelArray BigWheel7::wheel;
    inline BigWheel7::tpWheelMatrix BigWheel7::mods;
    inline std::array<BigWheel7::tpWheelArray*, BigWheel7::wheel_values> BigWheel7::factor;
}

// bits in a byte, allways EIGHT
constexpr unsigned byte_bits = 8; static_assert(byte_bits == 8);

// the batch is divided in chunks
// chuncks in one batch are processed sequentially / incrementally
// data for one chunk should fit in cache
constexpr unsigned chunk_bytes_lg2 = 18;   // 2^16 = 64k
constexpr unsigned chunk_bytes = 1 << chunk_bytes_lg2;
// try to have sizes as 2^n everytime
static_assert(chunk_bytes == 1 << chunk_bytes_lg2);
constexpr uint64_t chunk_bits = chunk_bytes * byte_bits;

// the whole interval is divided in segments
// each segment but the last one should fit N full chunks, with no rest
constexpr unsigned max_segment_bytes = 1 << 30; // 1 << 30 = 1GB