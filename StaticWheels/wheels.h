#pragma once

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

    constexpr auto range(unsigned max) { return std::views::iota(0u, max); }
    constexpr auto range(unsigned start, unsigned max) { return std::views::iota(start, max); }

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
            for (auto i : range(wheel_values))
                for (auto j : range(wheel_values))
                    mods[i][j] = wheel[i] * wheel[j] % wheel_span;
        }

        void GenFactors(void)
        {
            if (factors_initialized) return;
            cTimer tmr; tmr.Start();
            GenMods();
            unsigned values = 0;

            auto compute_line = [&](const unsigned i)
            {
                unsigned ni = wheel[i];
                for (auto j : range(wheel_values))
                {
                    unsigned nj = wheel[j];
                    for (auto k : range(wheel_values))
                    {
                        if (mods[i][k] == nj)
                        {
                            factor[i][j] = wheel[k];
                            //values++; 
                            break;
                        }
                    }
                }
            };
            std::for_each(std::execution::par_unseq, 
                range(wheel_values).begin(), range(wheel_values).end(), compute_line);

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
            tmr.Stop(true, "Factors initialized");
        }

        RootSieve() { GenFactors(); };
    };

//#ifdef TESTWHEELSINITS
    bool RootSieve<3>::factors_initialized = false;
    RootSieve<3>::tpWheelMatrix RootSieve<3>::factor;
    RootSieve<3>::tpWheelMatrix RootSieve<3>::mods;
    bool RootSieve<4>::factors_initialized = false;
    RootSieve<4>::tpWheelMatrix RootSieve<4>::factor;
    RootSieve<4>::tpWheelMatrix RootSieve<4>::mods;
    bool RootSieve<5>::factors_initialized = false;
    RootSieve<5>::tpWheelMatrix RootSieve<5>::factor;
    RootSieve<5>::tpWheelMatrix RootSieve<5>::mods;
    bool RootSieve<6>::factors_initialized = false;
    RootSieve<6>::tpWheelMatrix RootSieve<6>::factor;
    RootSieve<6>::tpWheelMatrix RootSieve<6>::mods;
//#endif
}
