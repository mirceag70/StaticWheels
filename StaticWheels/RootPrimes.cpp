#include "Helper.h"

inline uint8_t* sieve1; inline uint8_t* sieve5; inline uint8_t* sieve7; inline uint8_t* sieve11;

constexpr unsigned STEP_SIZE = 65'000'000;

//inline bool GetBit(tpPrime n, uint8_t sieve[])
//{
//	tpPrime q = (n / 12) / 8, r = (n / 12) % 8;
//	return (sieve[q] & BIT_MASK[r]);
//};
//inline void ResetBit(tpPrime n, uint8_t sieve[])
//{
//	tpPrime q = (n / 12) / 8, r = (n / 12) % 8;
//	sieve[q] &= BIT_RESET_MASK[r];
//};
//inline void FlipBit(tpPrime n, uint8_t sieve[])
//{
//	tpPrime q = (n / 12) / 8, r = (n / 12) % 8;
//	sieve[q] ^= BIT_MASK[r];
//};

inline void Pattern11(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 4));
	for (tpPrime x = 1, jmp = 0; x <= xmax; x += (1 + jmp), jmp = 1 - jmp)
	{
		//get in position
		tpPrime y, n0 = 4 * x * x;
		if (n0 < start)
		{
			const tpPrime yy = (tpPrime)ceil(sqrt(start - n0));
			y = 6 * (yy / 6) + 3;
			if (y < yy) y += 6;
		}
		else y = 3;
		//sieve
		for (; ; y += 6)
		{
			tpPrime n = n0 + (y * y);
			if (n > stop) break;
			FlipBit(n/12, sieve);
		}
	}
};
inline void Pattern12(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	//for x = 0 we can not start with y = 1
	uint8_t step = 2;
	const tpPrime yy = (tpPrime)ceil(sqrt(start));
	tpPrime y = 6 * (yy / 6) + 5;
	if (y < yy)
	{
		y += 2;
		if (y < yy) y += 4;
		else step = 4;
	}
	for (; ; y += step, step = 6 - step)
	{
		tpPrime n = y * y;
		if (n > stop) break;
		FlipBit(n/12, sieve);
	}
	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 4));
	for (tpPrime x = 3; x <= xmax; x += 3)
	{
		//get in position
		uint8_t step = 4;
		tpPrime y, n0 = 4 * x * x;
		if (n0 < start)
		{
			const tpPrime yy = (tpPrime)ceil(sqrt(start - n0));
			y = 6 * (yy / 6) + 1;
			if (y < yy)
			{
				y += 4;
				if (y < yy) y += 2;
				else step = 2;
			}
		}
		else y = 1;
		//sieve
		for (; ; y += step, step = 6 - step)
		{
			tpPrime n = (4 * x * x) + (y * y);
			if (n > stop) break;
			FlipBit(n/12, sieve);
		}
	}
};
inline void Pattern5(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 4));
	for (tpPrime x = 1, jmp = 0; x <= xmax; x += 1 + jmp, jmp = 1 - jmp)
	{
		//get in position
		uint8_t step = 4;
		tpPrime y, n0 = 4 * x * x;
		if (n0 < start)
		{
			const tpPrime yy = (tpPrime)ceil(sqrt(start - n0));
			y = 6 * (yy / 6) + 1;
			if (y < yy)
			{
				y += 4;
				if (y < yy) y += 2;
				else step = 2;
			}
		}
		else y = 1;
		//sieve
		for (; ; y += step, step = 6 - step)
		{
			tpPrime n = n0 + (y * y);
			if (n > stop) break;
			FlipBit(n/12, sieve);
		}
	}
};
inline void Pattern7(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 3));
	for (tpPrime x = 1; x <= xmax; x += 2)
	{
		//get in position
		uint8_t step = 2;
		tpPrime y, n0 = 3 * x * x;
		if (n0 < start)
		{
			const tpPrime yy = (tpPrime)ceil(sqrt(start - n0));
			y = 6 * (yy / 6) + 2;
			if (y < yy)
			{
				y += 2;
				if (y < yy) y += 4;
				else step = 4;
			}
		}
		else y = 2;
		//sieve			
		for (; ; y += step, step = 6 - step)
		{
			tpPrime n = n0 + (y * y);
			if (n > stop) break;
			FlipBit(n/12, sieve);
		}
	}
};
inline void Pattern111(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 2));
	tpPrime x = (tpPrime)ceil(sqrt(start / 3));
	if (x & 1) x++;
	for (; x <= xmax; x += 2)
	{
		//get in position
		uint8_t step = 4;
		tpPrime y, n0 = 3 * x * x;
		if (n0 > stop)
		{
			const tpPrime yy = (tpPrime)ceil(sqrt(n0 - stop));
			y = 6 * (yy / 6) + 1;
			if (y < yy)
			{
				y += 4;
				if (y < yy) y += 2;
				else step = 2;
			}
		}
		else y = 1;
		//sieve			
		for (; x > y; y += step, step = 6 - step)
		{
			tpPrime n = n0 - (y * y);
			if (n < start) break;
			FlipBit(n/12, sieve);
		}
	}
};
inline void Pattern112(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 2));
	tpPrime x = (tpPrime)ceil(sqrt(start / 3));
	if (not (x & 1)) x++;
	for (; x <= xmax; x += 2)
	{
		//get in position
		uint8_t step = 2;
		tpPrime y, n0 = 3 * x * x;
		if (n0 > stop)
		{
			const tpPrime yy = (tpPrime)ceil(sqrt(n0 - stop));
			y = 6 * (yy / 6) + 2;
			if (y < yy)
			{
				y += 2;
				if (y < yy) y += 4;
				else step = 4;
			}
		}
		else y = 2;
		//sieve			
		for (; x > y; y += step, step = 6 - step)
		{
			tpPrime n = n0 - (y * y);
			if (n < start) break;
			FlipBit(n/12, sieve);
		}
	}
};

template<typename FN1, typename FN2>
void Sieve2Patterns(FN1 fn1, FN2 fn2, const tpPrime limit, uint8_t sieve[])
{
	tpPrime iter_min = 0;
	if (limit > STEP_SIZE)
		for (; iter_min < (limit - STEP_SIZE); iter_min += STEP_SIZE)
		{
			fn1(iter_min, iter_min + STEP_SIZE, sieve);
			fn2(iter_min, iter_min + STEP_SIZE, sieve);
		}
	fn1(iter_min, limit, sieve); fn2(iter_min, limit, sieve);
};

template<typename FN>
void Sieve1Pattern(FN fn1, const tpPrime limit, uint8_t sieve[])
{
	tpPrime iter_min = 0;
	if (limit > STEP_SIZE)
		for (; iter_min < (limit - STEP_SIZE); iter_min += STEP_SIZE)
		{
			fn1(iter_min, iter_min + STEP_SIZE, sieve);
		}
	fn1(iter_min, limit, sieve);
};

void Sieve1(const tpPrime limit)
{
	Sieve2Patterns(Pattern11, Pattern12, limit, sieve1);
}
void Sieve5(const tpPrime limit)
{
	Sieve1Pattern(Pattern5, limit, sieve5);
}
void Sieve7(const tpPrime limit)
{
	Sieve1Pattern(Pattern7, limit, sieve7);
}
void Sieve11(const tpPrime limit)
{
	Sieve2Patterns(Pattern111, Pattern112, limit, sieve11);
}

void MarkSquares(uint8_t work_sieve[], const tpPrime limit, const unsigned rem)
{
	auto MarkOneSquare = [&](uint8_t sieve[], const tpPrime n)
	{
		if (GetBit(n/12, sieve))
			for (tpPrime i = rem * n * n; i < limit; i += 12 * n * n)
				ResetBit(i/12, work_sieve);
	};

	tpPrime kmax = (tpPrime)floor(sqrt(limit));

	MarkOneSquare(sieve5, 5);
	MarkOneSquare(sieve7, 7);
	MarkOneSquare(sieve11, 11);
	for (tpPrime k = 12; k < kmax; k += 12)
	{
		MarkOneSquare(sieve1, k + 1);
		MarkOneSquare(sieve5, k + 5);
		MarkOneSquare(sieve7, k + 7);
		MarkOneSquare(sieve11, k + 11);
	}
};

uint8_t svrp[4 * (ROOT_LIMIT / 96 + 2)];
uint8_t root_primes_gap[NO_ROOT_PRIMES+1];
constexpr unsigned NTC = 6; static_assert(NTC == 6, "NTC must be 6");
constexpr unsigned firstIdx7[NTC] = { 0, 125'793, 239'115, 348'510, 455'379, 560'685 };
constexpr unsigned gap7[NTC] = { 0, 30, 2, 36, 62, 30 };
constexpr unsigned firstIdx8[NTC] = { 0, 1'071'203, 2'050'938, 3'001'131, 3'933'307, 4'852'492 };
constexpr unsigned gap8[NTC] = { 0, 38, 36, 70, 60, 58 };
constexpr unsigned firstIdx9[NTC] = { 0, 9'327'000, 17'955'276, 26'355'864, 34'614'424, 42'770'642 };
constexpr unsigned gap9[NTC] = { 0, 24/2, 6/2, 62/2, 12/2, 28/2 };

const unsigned* firstIdx =  (ROOT_LIMIT_lg10 == 7) ? firstIdx7 : 
							(ROOT_LIMIT_lg10 == 8) ? firstIdx8 :
							(ROOT_LIMIT_lg10 == 9) ? firstIdx9 : NULL;
const unsigned* correction_gap = (ROOT_LIMIT_lg10 == 7) ? gap7 :
								 (ROOT_LIMIT_lg10 == 8) ? gap8 :
								 (ROOT_LIMIT_lg10 == 9) ? gap9 : NULL;

inline thread_local unsigned idx_last_prime_t = 0;
inline thread_local unsigned last_prime_t = 0;

inline void AddRootPrime(unsigned prime)
{
#if LIMIT_lg10 < 9
	//assert(prime - last_prime_t < 256);
	root_primes_gap[idx_last_prime_t++] = (uint8_t)(prime - last_prime_t);
#else
	//assert(prime - last_prime_t < 512);
	root_primes_gap[idx_last_prime_t++] = (uint8_t)((prime - last_prime_t) / 2);
#endif
	last_prime_t = prime;
}
void DoOneChunk(const unsigned svlen, const unsigned nt, const unsigned i, unsigned* numprimes)
{
	const unsigned chunksz = svlen / nt;
	const unsigned start = i * chunksz;
	const unsigned stop = (i == (nt - 1)) ? svlen : start + chunksz;

	idx_last_prime_t = firstIdx[i];

	auto CountPrime = [&](const unsigned n) { (*numprimes)++; AddRootPrime(n); };
	for (unsigned k = start; k < stop; k++)
	{
		uint8_t sv1 = sieve1[k], sv5 = sieve5[k], sv7 = sieve7[k], sv11 = sieve11[k];

		unsigned n = 12 * (8 * k);
		if (sv1 & BIT_MASK[0]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[0]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[0]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[0]) CountPrime(n + 11);
		n = 12 * (8 * k + 1);
		if (sv1 & BIT_MASK[1]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[1]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[1]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[1]) CountPrime(n + 11);
		n = 12 * (8 * k + 2);
		if (sv1 & BIT_MASK[2]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[2]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[2]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[2]) CountPrime(n + 11);
		n = 12 * (8 * k + 3);
		if (sv1 & BIT_MASK[3]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[3]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[3]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[3]) CountPrime(n + 11);
		n = 12 * (8 * k + 4);
		if (sv1 & BIT_MASK[4]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[4]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[4]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[4]) CountPrime(n + 11);
		n = 12 * (8 * k + 5);
		if (sv1 & BIT_MASK[5]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[5]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[5]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[5]) CountPrime(n + 11);
		n = 12 * (8 * k + 6);
		if (sv1 & BIT_MASK[6]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[6]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[6]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[6]) CountPrime(n + 11);
		n = 12 * (8 * k + 7);
		if (sv1 & BIT_MASK[7]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[7]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[7]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[7]) CountPrime(n + 11);
	}
}

//root primes generator
void SoA_LP_gen_root_primes(void)
{
	const unsigned svlen = (ROOT_LIMIT / 96) + 2;
	sieve1 = svrp; sieve5 = sieve1 + svlen;
	sieve7 = sieve5 + svlen; sieve11 = sieve7 + svlen;

	// 2 and 3 are not generated by the alghorithm
	for (unsigned i = 0; i < 4 * svlen; i++) svrp[i] = false;

	// sieve quadratics
	{
		std::thread t1 = std::thread(&Sieve1, ROOT_LIMIT);
		std::thread t5 = std::thread(&Sieve5, ROOT_LIMIT);
		std::thread t7 = std::thread(&Sieve7, ROOT_LIMIT);
		Sieve11(ROOT_LIMIT);
		t7.join(); t5.join(); t1.join();
	}

	// sieve multiple of squares
	{
		std::thread t1 = std::thread(&MarkSquares, sieve1, ROOT_LIMIT, 1);
		std::thread t5 = std::thread(&MarkSquares, sieve5, ROOT_LIMIT, 5);
		std::thread t7 = std::thread(&MarkSquares, sieve7, ROOT_LIMIT, 7);
		MarkSquares(sieve11, ROOT_LIMIT, 11);
		t7.join(); t5.join(); t1.join();
	}

	// get primes
	unsigned numprimes = 0;
	{
		std::vector<std::thread> tc;
		const unsigned nt = NTC; unsigned np[nt]{ 0 };
		for (auto i : Range<nt>())
			tc.push_back(std::thread(&DoOneChunk, svlen, nt, i, np + i));
		for (auto& t : tc) t.join();
		for (auto i : Range<nt>()) numprimes += np[i];
		assert(numprimes == NO_ROOT_PRIMES);
	}

	//fix values
	for (unsigned i = 1; i <= 5; i++) root_primes_gap[firstIdx[i]] = correction_gap[i];

#define	TEST_ROOT_PRIMES false
	// test root primes
#if TEST_ROOT_PRIMES == true
	AddPrime(2); AddPrime(3); AddPrime(5);
	uint64_t p = 5;
	for (unsigned i = 1; i < NO_ROOT_PRIMES; i++)
	{
#if LIMIT_lg10 < 9
		p += root_primes_gap[i];
#else
		p += 2 * root_primes_gap[i];
#endif
		AddPrime(p);
	}
	std::cout << "Root Primes initialization success for N = " << LIMIT_lg10;
	exit(15);
#endif

	root_primes_gap[NO_ROOT_PRIMES] = 255;
}
