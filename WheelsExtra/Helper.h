#pragma once

#include <iostream>
#include <iomanip>
#include <cassert>
#include <array>
#include <thread>
#include <numeric>
#include <fstream>
#include <sstream>
#include <functional>

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
		if(not msg.empty()) std::cout << "[" << msg <<"] ";
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

	void Stop(bool const print_time = false)
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
			OutputTime(time_accumulated, "StopTime");
	}

	long LapTime(bool print_time = false) const
	{
		assert(timer_state == InternalStates::Counting);
		
		// get current clock
		clock_t tm = clock();
		long laptime = TimeDiff2Miliseconds(tm - mark_lap);
		//prepare for next lap
		mark_lap = tm;
		
		if (print_time)
			OutputTime(laptime, "LapTime");
		
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

template<std::size_t range_length, typename INTEGRAL = int>
consteval auto Range(const int offset = 0)
{
	static_assert(std::is_integral<INTEGRAL>::value, 
				"Integral type required for Range.");

	std::array<INTEGRAL, range_length> arr{};
	for (int i = 0; i < range_length; i++)
		arr[i] = i + offset;
	return arr;
};

inline const void nln(bool delimiter = false)
{
	std::cout << std::endl;
	if (delimiter)
		std::cout << "\t--------------------------" << std::endl;
}

template<class NUMBER>
constexpr double Average(std::vector<NUMBER> vec)
{
	static_assert(std::is_arithmetic<NUMBER>::value,
		"Numeric type required for Average.");

	return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

// scalar data type used for prime numbers
//go with 32b for experiments
#define USE_64_BITS_PRIMES
//we need 64b for larger numbers
#ifndef USE_64_BITS_PRIMES
typedef		uint32_t		tpPrime;
#else
typedef		uint64_t		tpPrime;
#endif // !USE_64_BITS

typedef		uint32_t		tpRoot;


class cWriter
{
private:
	static constexpr auto no_per_row = 10;
	static constexpr auto no_per_file = 100'000;
	static_assert(no_per_file % no_per_row == 0, "Must have full lines per output file");

	static constexpr auto counter_places = 5;	//counter - suffix for filename
	static constexpr auto file_extension = ".primes";

	std::string file_name;
	std::ofstream file_with_primes;
	
	unsigned idx_number = 0; 
	unsigned idx_file = 1;

	void OpenNewFile(void)
	{
		std::ostringstream oss(file_name, std::ostringstream::ate);
		//format file name
		oss << std::setw(counter_places) << std::setfill('0') << idx_file;
		oss << file_extension;

		assert(!file_with_primes.is_open());
		file_with_primes.open(oss.str(), std::ifstream::out | std::fstream::trunc);
		assert(file_with_primes.is_open());

		idx_file++;
	}

	void CloseFile(void)
	{
		file_with_primes.flush();
		file_with_primes.close();
		assert(!file_with_primes.is_open());
	}

public:
	void WritePrime(tpPrime n)
	{
		idx_number++;

		//new file if necessary
		if (idx_number > no_per_file)
		{
			CloseFile();
			OpenNewFile();
			idx_number = 1;
		}

		//write it
		file_with_primes << std::setw(15) << n;
		if (idx_number % no_per_row == 0)
			file_with_primes << std::endl;
	}

	cWriter(const char* fName) :
		file_name(fName)
	{
		OpenNewFile();
	}

	~cWriter()
	{
		CloseFile();
	}
};

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

/*
// C++14 compile time square root using binary search
template <typename T>
constexpr T sqrt_helper(T x, T lo, T hi)
{
	if (lo == hi)
		return lo;

	const T mid = (lo + hi + 1) / 2;

	if (x / mid < mid)
		return sqrt_helper<T>(x, lo, mid - 1);
	else
		return sqrt_helper(x, mid, hi);
}

template <typename INTEGRAL>
constexpr INTEGRAL ConstSqrtT(INTEGRAL x)
{
	static_assert(std::is_integral<INTEGRAL>::value,
		"Integral type required for ConstSqrt.");

	return sqrt_helper<INTEGRAL>(x, 0, x / 2 + 1);
}
*/

template<class NUMBER>
constexpr NUMBER ConstAbs(NUMBER a)
{
	static_assert(std::is_arithmetic<NUMBER>::value,
		"Numeric type required for ConstAbs.");

	return (a >= 0) ? a : (-a);
}

template<class NUMBER>
constexpr long long ConstRound(NUMBER a)
{
	static_assert(std::is_arithmetic<NUMBER>::value,
		"Numeric type required for ConstRound.");

	//if (a == 0)
		//return 0;

	long long rnd = static_cast<long long>(a);
	double delta = static_cast<double>(ConstAbs(a - rnd));

	if (delta <= 0.50)
		return rnd;
	else
		return rnd > 0 ? rnd + 1 : rnd - 1;
}

static_assert(0 == ConstRound(0));
static_assert(5 == ConstRound(5.35));
static_assert(-5 == ConstRound(-5.35));
static_assert(6 == ConstRound(5.75));
static_assert(-6 == ConstRound(-5.65));

// square root of an integral number using Newtons method
template <typename INTEGRAL>
consteval INTEGRAL ConstSqrt(INTEGRAL n)
{
	static_assert(std::is_integral<INTEGRAL>::value,
		"Integral type required for ConstSqrt.");

	double x = n / 10.0;			// current value
	double gss = (x + (n / x)) / 2;	// current guess 

	while(ConstAbs(gss - x) > 0.25)
	// we need an integer, so our tolerance is quite high
	{
		x = gss;
		gss = (x + (n / x)) / 2;	// better guess
	}

	return static_cast<INTEGRAL>(ConstRound(gss));
}

static_assert(5 == ConstSqrt(25));
static_assert(27345 == ConstSqrt(27345 * 27345));
static_assert(1'000'000 == ConstSqrt(1'000'000'000'000));

#define	USE_ARITHMETIC_LOGIC

class cTriState
{
public:
	enum class State : char { Negative = -1, Zero = 0, Positive = 1};

	State GetState(void) { return inner_state; };
	void SetState(State newState) { inner_state = newState; };

	void IncState(void) 
	{
#ifdef USE_ARITHMETIC_LOGIC
		inner_state = static_cast<State>(static_cast<char>(inner_state) + 1 - 3 * (inner_state == State::Positive));
#else
		switch (inner_state)
		{
		case State::Negative:
			inner_state = State::Zero;
			break;
		case State::Zero:
			inner_state = State::Positive;
			break;
		case State::Positive:
			inner_state = State::Negative;
			break;
		}
		return inner_state; 
#endif // USE_ARITHMETIC_LOGIC
	};
	
	void DecState(void) 
	{
#ifdef USE_ARITHMETIC_LOGIC
		inner_state = static_cast<State>(static_cast<char>(inner_state) - 1 + 3 * (inner_state == State::Negative));
#else
		switch (inner_state)
		{
		case State::Negative:
			inner_state = State::Positive;
			break;
		case State::Zero:
			inner_state = State::Positive;
			break;
		case State::Positive:
			inner_state = State::Zero;
			break;
		}
		return inner_state;
#endif // USE_ARITHMETIC_LOGIC
	};

	friend bool operator>(cTriState& tri, State state) { return tri.inner_state > state; };
	friend bool operator==(cTriState& tri, State state) { return tri.inner_state == state; };
	
	void operator++(int) { IncState(); }
	void operator=(State state) { SetState(state); }

	cTriState(State newState) { inner_state = newState; };
	cTriState(void) {};

private:
	State inner_state = State::Zero;
};

consteval cTriState::State operator "" _tri(const char c)
{
	switch (c)
	{
	case 'N': return cTriState::State::Negative;
	case 'Z': return cTriState::State::Zero;
	case 'P': return cTriState::State::Positive;
	}
}

static_assert(cTriState::State::Negative < cTriState::State::Zero);

constexpr unsigned firstPrimes[]{
	  2,    3,    5,    7,   11,   13,   17,   19,   23,  29,
	 31,   37,   41,   43,   47,   53,   59,   61 ,  67,  71,
	 73,   79,   83,   89,   97,  101,  103,  107,  109,  113,  
	127,  131,  137,  139,  149,  151,  157,  163,  167,  173,  
	179,  181,  191,  193,  197,  199,  211,  223,  227,  229,  
	233,  239,  241,  251,  257,  263,  269,  271,  277,  281,  
	283,  293,  307,  311,  313,  317,  331,  337,  347,  349,  
	353,  359,  367,  373,  379,  383,  389,  397,  401,  409,  
	419,  421,  431,  433,  439,  443,  449,  457,  461,  463,  
	467,  479,  487,  491,  499,  503,  509,  521,  523,  541,  
	547,  557,  563,  569,  571,  577,  587,  593,  599,  601
};

template<class T>
void DeleteVectorElements(std::vector<T>& v)
{
	//for (; not v.empty(); v.pop_back()) delete v.back();
}

consteval auto BuilMaskArray(void)
{
	std::array<unsigned char, 8> arr{};
	for (unsigned i = 0; i < 8; i++)
		arr[i] = (1 << (i));
	return arr;
}

inline auto v_msk = BuilMaskArray();

#define BIT_ARRAY
#ifdef BIT_ARRAY

const unsigned MULT_FACTOR = 12*8;

template<typename T1, typename T2>
constexpr auto FLIP12(T1 vect, T2 pos)
{
	T2 q = (pos / 12) / 8;
	T2 r = (pos / 12) % 8;
	vect[q] ^= v_msk[r];// ^ vect[q];
}

template<typename T1, typename T2>
constexpr auto GET12(T1 vect, const T2 pos)
{
	T2 q = (pos / 12) / 8;
	T2 r = (pos / 12) % 8;
	return vect[q] & v_msk[r];
}

template<typename T1, typename T2>
constexpr auto RESET12(T1 vect, const T2 pos)
{
	T2 q = (pos / 12) / 8;
	T2 r = (pos / 12) % 8;
	vect[q] &= ~v_msk[r];
}

#else

const unsigned MULT_FACTOR = 12;

template<typename T1>
constexpr auto FLIP12(T1 vect, const size_t pos) { 
	vect[(pos) / 12] ^= true; 
}

template<typename T1>
constexpr auto GET12(const T1 vect, const size_t pos) { return vect[(pos) / 12]; }

template<typename T1>
constexpr auto RESET12(T1 vect, const size_t pos) { 
	vect[pos / 12] = false; 
}

#endif // BIT_ARRAY

inline uint64_t PI_Nmax(uint64_t Nmax)
{
	double logNmax = log(Nmax);
	return (uint64_t)(Nmax / logNmax * (1 + 1.2762 / logNmax));
}

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

	last_primes[idx_last_primes++] = prime;
}

inline thread_local uint8_t idx_last_primes_t = 0;
inline thread_local std::array<uint64_t, 256> last_primes_t;
inline void AddPrimeThreaded(uint64_t prime)
{
	last_primes_t[idx_last_primes_t++] = prime;
}

typedef std::array<uint8_t, 256> tp256BitV;
consteval tp256BitV generate_bit1_count_table()
{
	tp256BitV v = { 0 };
	//bit 1 count
	for (int i = 1; i < v.size(); i++) v[i] = v[i / 2] + (i & 1);
	return v;
}
consteval tp256BitV generate_bit0_count_table()
{
	tp256BitV v = { 0 };
	//bit 1 count
	for (int i = 1; i < v.size(); i++) v[i] = v[i / 2] + (i & 1);
	//bit 0 count 
	for (int i = 0; i < v.size(); i++) v[i] = 8 - v[i];
	return v;
}
constexpr tp256BitV bit1_count_table = generate_bit1_count_table();
constexpr tp256BitV bit0_count_table = generate_bit0_count_table();