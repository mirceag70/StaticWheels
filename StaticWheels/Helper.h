#pragma once

#include <iostream>
#include <iomanip>
#include <queue>
#include <cassert>
#include <array>
#include <thread>
#include <future>
#include <numeric>
#include <execution>
#include <fstream>
#include <chrono>
#include <limits>
#include <functional>
#include <typeinfo>
#include <ranges>

class cTimer
{
private:

	enum class InternalStates { Initial, Iddle, Counting, Pause, Stop };

	InternalStates timer_state = InternalStates::Initial;

	long time_accumulated = 0;
	clock_t mark_start = 0;
	mutable clock_t mark_lap = 0;

	constexpr long TimeDiff2Miliseconds(clock_t const tmDiff) const
	{
		return tmDiff * (1'000 / CLOCKS_PER_SEC);
	}

	void OutputTime(long const tm, std::string_view const msg = "") const
	{
		std::cout << " ( ";
		if (not msg.empty()) std::cout << "[" << msg << "] ";
		std::cout << std::fixed << std::setprecision(0) << tm;
		std::cout << " ms ) ";
	}

public:

	cTimer(void) : timer_state(InternalStates::Iddle) {};

	void Start(void)
	{
		assert(timer_state == InternalStates::Iddle);
		timer_state = InternalStates::Counting;

		// mark down current clock
		mark_lap = mark_start = clock();	// use portable clock()
	}

	void Stop(bool const print_time = false, const char* msg = NULL)
	{
		switch (timer_state)
		{
		case InternalStates::Counting:
			time_accumulated += TimeDiff2Miliseconds(clock() - mark_start);
			mark_start = 0;
			[[fallthrough]];
		case InternalStates::Pause:
			timer_state = InternalStates::Stop;
			break;
		default:
			assert(false);
		}

		if (print_time)
			OutputTime(time_accumulated, msg ? msg : "StopTime");
	}

	long LapTime(const bool print_time = false, const char* msg = NULL) const
	{
		assert(timer_state == InternalStates::Counting);

		// get current clock
		clock_t tm = clock();
		long laptime = TimeDiff2Miliseconds(tm - mark_lap);
		//prepare for next lap
		mark_lap = tm;

		if (print_time)
			OutputTime(laptime, msg ? msg : "LapTime");

		return laptime;
	}

	long GetTime(bool print_time = false) const
	{
		clock_t tm = 0;

		switch (timer_state)
		{
		case InternalStates::Counting:
			tm = TimeDiff2Miliseconds(tm - mark_start);
			[[fallthrough]];
		case InternalStates::Pause:
		case InternalStates::Stop:
			tm += time_accumulated;
			break;
		default:
			assert(false);
		}

		if (print_time)
			OutputTime(tm, "Time");

		return tm;
	}

	long Pause(bool print_time = false)
	{
		assert(timer_state == InternalStates::Counting);
		timer_state = InternalStates::Pause;

		time_accumulated += TimeDiff2Miliseconds(clock() - mark_start);
		mark_start = 0;

		if (print_time)
			OutputTime(time_accumulated, "PauseTime");

		return time_accumulated;
	}

	void Resume(void)
	{
		assert(timer_state == InternalStates::Pause);
		timer_state = InternalStates::Counting;

		mark_lap = mark_start = clock();
	}

	void Reset(void)
	{
		assert(timer_state == InternalStates::Stop);
		timer_state = InternalStates::Iddle;

		time_accumulated = mark_lap = mark_start = 0;
	}
};

constexpr std::size_t default_iterations = 5;

inline const void nln(bool delimiter = false)
{
	std::cout << std::endl;
	if (delimiter)
		std::cout << "\t--------------------------" << std::endl;
}

template<std::size_t range_length, typename INTEGRAL = int>
consteval auto Range(const int offset = 0)
{
	static_assert(std::is_integral<INTEGRAL>::value,
		"Integral type required for Range.");

	std::array<INTEGRAL, range_length> arr{};
	for (int i = 0; i < range_length; i++)
		arr[i] = i + offset;
	return arr;
}

template<class NUMBER>
constexpr double Average(std::vector<NUMBER> vec)
{
	static_assert(std::is_arithmetic<NUMBER>::value,
		"Numeric type required for Average.");

	return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

inline uint64_t PI_Nmax(uint64_t Nmax)
{
	double logNmax = log(Nmax);
	return (uint64_t)(Nmax / logNmax * (1 + 1.2762 / logNmax));
}

// scalar data type used for prime numbers
//go with 32b for small experiments
#define USE_64_BITS_PRIMES
//we need 64b for larger numbers		
#ifndef USE_64_BITS_PRIMES
typedef		uint32_t		tpPrime;
#else
typedef		uint64_t		tpPrime;
#endif // !USE_64_BITS

class cChecker
{
	std::ifstream file_with_primes;
	tpPrime max_value = 0;

public:
	tpPrime check_next_prime(tpPrime nGenerated)
	{
		if (nGenerated == 2)
			BackToZero();

		tpPrime nRead;
		file_with_primes >> nRead;
		if (nRead <= max_value)
			return nRead == nGenerated;
		else
			return true;
	}

	void BackToZero(void)
	{
		file_with_primes.seekg(0, file_with_primes.beg);
	}

	cChecker(std::string file_name = "C:/50MilPrimes.txt", tpPrime maxPrime = 982'451'653) :
		max_value(maxPrime)
	{
		file_with_primes.open(file_name, std::ifstream::in);
		assert(file_with_primes.is_open());
	}

	~cChecker() { file_with_primes.close(); }
};

inline uint8_t idx_last_primes = 0;
inline std::array<uint64_t, 256> last_primes;
//inline uint32_t maxgap, biggap;
//inline uint64_t lastprime;
inline void AddPrime(uint64_t prime)
{ 
#ifdef _DEBUG
	static cChecker ckr;
	assert(ckr.check_next_prime(prime));
#endif // _DEBUG
	
	//if (prime == 2) maxgap = biggap = 0;
	//else
	//if ((prime - lastprime) > maxgap)
	//{
	//	maxgap = (uint32_t)(prime - lastprime);
	//	if (maxgap > 0xFF) biggap++;
	//}
	//lastprime = prime;

	last_primes[idx_last_primes++] = prime; 
}

inline thread_local uint8_t idx_last_primes_t = 0;
inline thread_local std::array<uint64_t, 256> last_primes_t;
inline void AddPrimeThreaded(uint64_t prime)
{
	last_primes_t[idx_last_primes_t++] = prime;
}

#ifndef ITERATIONS
#define ITERATIONS	5
#endif 

constexpr unsigned firstPrimes[] = {
		 2,         3,         5,         7,        11,        13,        17,        19,
		23,        29,        31,        37,        41,        43,        47,        53,
		59,        61,        67,        71,        73,        79,        83,        89,
		97,       101,       103,       107,       109,       113,       127,       131,
	   137,       139,       149,       151,       157,       163,       167,       173,
	   179,       181,       191,       193,       197,       199,       211,       223,
	   227,       229,       233,       239,       241,       251,       257,       263,
	   269,       271,       277,       281,       283,       293,       307,       311,
	   313,       317,       331,       337,       347,       349,       353,       359,
	   367,       373,       379,       383,       389,       397,       401,       409,
	   419,       421,       431,       433,       439,       443,       449,       457,
	   461,       463,       467,       479,       487,       491,       499,       503,
	   509,       521,       523,       541,       547,       557,       563,       569,
	   571,       577,       587,       593,       599,       601,       607,       613,
	   617,       619,       631,       641,       643,       647,       653,       659,
	   661,       673,       677,       683,       691,       701,       709,       719,
	   727,       733,       739,       743,       751,       757,       761,       769,
	   773,       787,       797,       809,       811,       821,       823,       827,
	   829,       839,       853,       857,       859,       863,       877,       881,
	   883,       887,       907,       911,       919,       929,       937,       941,
	   947,       953,       967,       971,       977,       983,       991,       997 };

template<typename INTEGRAL1, typename INTEGRAL2, typename INTEGRAL3,
		std::size_t size1, std::size_t size2 = 0, std::size_t size3 = 0>
void Try_Sieve(const uint64_t LIMIT, std::string message,
	std::function<uint64_t(uint64_t, INTEGRAL1[], INTEGRAL2[], INTEGRAL3[])> Sieve, bool show_all = true)
{
	static_assert(std::is_integral<INTEGRAL1>::value, "Integral type required for INTEGRAL1.");
	static_assert(std::is_integral<INTEGRAL2>::value, "Integral type required for INTEGRAL2.");
	static_assert(std::is_integral<INTEGRAL3>::value, "Integral type required for INTEGRAL3.");

	INTEGRAL1* v1 = NULL; INTEGRAL2* v2 = NULL; INTEGRAL3* v3 = NULL;

	if (size1 > 0) v1 = new INTEGRAL1[size1];
	if (size2 > 0) v2 = new INTEGRAL2[size2];
	if (size3 > 0) v3 = new INTEGRAL3[size3];

	nln(true);
	cTimer tmr;
	std::vector<decltype(tmr.LapTime())> times;

	tmr.Start();
	std::cout << " - " << message << " - ";
	for (auto i : Range<ITERATIONS>())
	{

		auto numPrimes = Sieve(LIMIT, v1, v2, v3);

		if (i == 0 or show_all)
		{
			nln();
			std::cout << numPrimes << " primes up to " << LIMIT;
			times.push_back(tmr.LapTime(true));
		}
		else
		{
			times.push_back(tmr.LapTime(false));
		}
	}
	nln(true);
	tmr.Stop();
	std::cout << "Average compute time: " << Average(times);

	if (v1) delete[] v1;
	if (v2) delete[] v2;
	if (v3) delete[] v3;
}

constexpr uint8_t BIT_MASK[] = { 1, 2, 4, 8, 16, 32, 64, 128 };
constexpr uint8_t BIT_RESET_MASK[] = { ~1, ~2, ~4, ~8, ~16, ~32, ~64, ~128 };

constexpr void SetBit(const uint64_t bitidx, uint8_t sv[])
{
	uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
	sv[idx] |= BIT_MASK[bit];
};
constexpr void ResetBit(uint64_t bitidx, uint8_t sv[])
{
	uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
	sv[idx] &= BIT_RESET_MASK[bit];
};
constexpr auto GetBit(uint64_t bitidx, uint8_t sv[])
{
	uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
	return (sv[idx] & BIT_MASK[bit]);
};
constexpr void FlipBit(uint64_t bitidx, uint8_t sv[])
{
	uint64_t idx = bitidx / 8; uint8_t bit = (uint8_t)(bitidx % 8);
	sv[idx] ^= BIT_MASK[bit];
};

inline uint64_t idx2no(uint64_t idx) { return 3 * idx + 5 - (idx & 1); };
inline uint64_t no2idx(uint64_t no) { return no / 3 - 1; };

// use 7 up to 10^14, 8 up to 10^16, 9 up tp 10^18 or 10 up to 2^64
// for 9 and 10 the gaps are halfed
#define LIMIT_lg10 7
constexpr unsigned ROOT_LIMIT_lg10 = LIMIT_lg10; 
static_assert( ROOT_LIMIT_lg10 >= 7 or ROOT_LIMIT_lg10 <= 9, "ROOT_LIMIT must be 7, 8 or 9" );
//10 is not yet corrected

constexpr unsigned PI[11] = // minus 2, for primes 2 and 3 (which are not computed)
	{ 0,0,0,0,0,0,    78'498,    664'577,   5'761'453,     50'847'532, 203'280'219 };
constexpr unsigned ROOT_LIMITS[11] = 
	{ 0,0,0,0,0,0, 1'000'000, 10'000'000, 100'000'000,  1'000'000'000, 0xFFFF'FFFF };

constexpr unsigned NO_ROOT_PRIMES = PI[ ROOT_LIMIT_lg10 ];
constexpr unsigned ROOT_LIMIT = ROOT_LIMITS[ROOT_LIMIT_lg10 ];

// for code supposed to be unreachable
#ifdef DEBUG
# define NEVERHERE   assert(0)
#else
# define NEVERHERE   __assume(0)
#endif

typedef std::array<uint8_t, 256> tp256BitV;
consteval tp256BitV generate_bit0_count_table()
{
	tp256BitV v = { 0 };
	//bit 1 count
	for (int i = 1; i < v.size(); i++) v[i] = v[i / 2] + (i & 1);
	//bit 0 count 
	for (int i = 0; i < v.size(); i++) v[i] = 8 - v[i];
	return v;
}
constexpr tp256BitV bit0_count_table = generate_bit0_count_table();