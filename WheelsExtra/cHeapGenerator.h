#pragma once

#include "Helper.h"

#include <thread>
#include <future>

/*
struct cPrime
{
public:
	// CurrentValue first, because it is the most accesses value
	tpPrime current_value;		// current value in mobile sieve
	//tpPrime prime_value;		// value of the prime number
	//tpRoot root_idx;			// current index in root sieve
	static constexpr unsigned bits_idx = 26;
	static constexpr unsigned mask_idx = (1 << bits_idx) - 1;
	tpPrime pck_val_idx;			//packed vaue and root idx
	static constexpr unsigned bits_val = 8 * sizeof(pck_val_idx) - bits_idx;
	//[val:38bits][idx:26bits]

	inline tpPrime GetVal(void)
	{
		return pck_val_idx >> bits_idx;
	}
	inline tpRoot GetIdx(void)
	{
		return pck_val_idx & mask_idx;
	}
	inline void SetValIdx(tpPrime val, tpRoot idx)
	{
		pck_val_idx = (val << bits_idx) + idx;
	}
	inline void IncIdx(void)
	{
		pck_val_idx++;
	}
	inline void SetIdx2Zero(void)
	{
		pck_val_idx &= ~mask_idx;
	}
};
*/

struct cPrimeVal
{
public:
	tpPrime prime_value;		// value of the prime number
	tpPrime root_idx;			// current index in root sieve
};

struct cPrimeHep
{
public:
	tpPrime current_value;		// current value in mobile sieve
	cPrimeVal* prime;			// handle for the prime
};

class cRootSieve
{
protected:
	static constexpr unsigned root_N = 9;
	static constexpr tpRoot root_sieve_size = 36'495'360;
	static constexpr tpRoot root_primorial = 223'092'870;
	typedef char rootSieveData;
	//typedef tpRoot rootSieveData;

	rootSieveData* root_sieve = nullptr;
	tpRoot* root_sieve_vals = nullptr;

	void InitializeRootSieve(void);
	void DeinitializeRootSieve(void) 
	{ 
		if (root_sieve) delete[] root_sieve; 
		if (root_sieve_vals) delete[] root_sieve_vals;
	}
};

class cHeapSieve : protected cRootSieve
{
protected:
	static constexpr unsigned row_size = 2 << 20; // 1'000'000;
	//each array in the vector has this fixed size
	//MUST be an EVEN number to assure idx_i for the 
	//last lift w/o sibling is ODD, for simpler tests
	static_assert(row_size % 2 == 0);

	static constexpr auto MAX_J_VECTOR_SIZE = 1'000;

#define USE_HEP_ROOT

#ifdef USE_HEP_ROOT
	cPrimeHep hep_root{};
	cPrimeVal val_root{};
#endif
	cPrimeHep* hep_primes[MAX_J_VECTOR_SIZE];
	cPrimeVal* val_primes[MAX_J_VECTOR_SIZE];

	inline void InitializeNode(cPrimeVal* pNV, cPrimeHep* pNH, tpPrime val);

	int idx_i = -1; //index of the last element in the last row
	int idx_j = -1; //index of the last row in vector
	//last element in the vector is at idx [i] in array [j]
	//having this address: vec_primes[idx_j]->data() + idx_i

	inline void AddNewRow(void);
	inline void AddNewNode(tpPrime val);

	void InitializeHeapSieve(void);
	void DeinitializeHeapSieve(void);
};

class cHeapGenerator : protected cHeapSieve
{
private:
	static constexpr unsigned no_of_children = 4;
	static constexpr unsigned manually_inserted_primes = 1 + no_of_children;

	cChecker checker;
	cWriter writer;

	tpPrime candidate = 1;
	tpPrime min_val_now;
	tpPrime limit = 0;
	tpPrime sqrtlimit = 0;

	std::promise<void> exitSignal;
	std::future<void> exitFuture;

	void Writer(std::future<void> futureObj);
	void WritePrimes(int offset = 0);
	void PrintStatus(void);

	inline void InsertNewPrime(void);
	inline void AdvanceSieve(void);
	inline void Process(void);
	inline void Heapify4(cPrimeHep OldRoot);
	inline void Heapify2(cPrimeHep OldRoot);
	inline void AdvanceIdxs4(tpPrime& idxi_child, tpPrime& idxj_child, const tpPrime idx_smlst);

	//case counters
	//tpPrime cnt_adv = 0;
	//tpPrime cnt_eq = 0;
	//tpPrime cnt_rsd = 0;
	//tpPrime cnt_swtchright = 0;

	tpPrime nr_hep_primes = 0;

public:
	cHeapGenerator(void);
	~cHeapGenerator(void);

	tpPrime no_gen_primes = 0;

	void Generator(tpPrime lmt);

	tpPrime GeneratePrimes(tpPrime lmt);
};
