#include "cHeapGenerator.h"

void cRootSieve::InitializeRootSieve(void)
{
    std::cout << "\nInitialize Root Sieve for N = 9";

    tpRoot i, idx = 0;  //counters

    cTimer tmr;
    tmr.Start();

    root_sieve_vals = new tpRoot[root_sieve_size];
    root_sieve = new rootSieveData[root_sieve_size];

    cTriState tri('P'_tri);
    for (i = 1; i < root_primorial; tri++, i += 2) //step over even numbers
        if (tri > 'N'_tri) // ignore multiple of 3
        {
            if ((i % 5) and (i % 7) and (i % 11) and (i % 13)
                and (i % 17) and (i % 19) and (i % 23))
                root_sieve_vals[idx++] = i;
        }

    assert(root_sieve_size == idx);

    // transform absolute values in step increments
    for (i = 1; i < root_sieve_size; i++)
#pragma warning(suppress:6001) // Using uninitialized memory '*root_sieve_tmp'
        root_sieve[i - 1] = root_sieve_vals[i] - root_sieve_vals[i - 1];
    root_sieve[root_sieve_size - 1] = root_primorial + 1 - root_sieve_vals[root_sieve_size - 1];

    tmr.Stop(true); nln();
}

void cHeapSieve::AddNewRow(void)
{
    idx_j++;
    assert(idx_j < MAX_J_VECTOR_SIZE);
    val_primes[idx_j] = (cPrimeVal*)_aligned_malloc(row_size * sizeof(cPrimeVal), 16);
    if(hep_primes[idx_j] = (cPrimeHep*)_aligned_malloc(row_size * sizeof(cPrimeHep), 64))
        std::memset(hep_primes[idx_j], 0xFF, row_size * sizeof(cPrimeHep));
    //std::cout << "new row;";
}

void cHeapSieve::DeinitializeHeapSieve(void)
{
#ifdef USE_STL_VECTOR
    //do nothing, vector auto-deletes
#else
    while(idx_j >= 0)
    {
        _aligned_free(hep_primes[idx_j]);
        _aligned_free(val_primes[idx_j]);
        idx_j--;
    }
#endif
}

void cHeapSieve::InitializeHeapSieve(void)
{
#ifdef USE_STL_VECTOR
    //vec_primes.resize(MAX_J_VECTOR_SIZE);
    hep_primes.resize(MAX_J_VECTOR_SIZE);
#endif
    AddNewRow();
}

void cHeapSieve::InitializeNode(cPrimeVal* pNV, cPrimeHep* pNH, tpPrime val)
{
    //set the new value
    pNV->prime_value = val;
    //advance one step
    pNH->prime = pNV;
    //get to sqr value
    
    //set initial value
    pNH->current_value = pNV->prime_value * pNV->prime_value;

    tpPrime mval = pNV->prime_value % root_primorial;
    //tpRoot* pix = std::lower_bound(root_sieve_vals, root_sieve_vals + root_sieve_size, mval);
    //assert(*pix == mval);
    //tpPrime ix = pix - root_sieve_vals;

    static tpPrime ixo = 1;
    while (root_sieve_vals[ixo] != mval)
        if (++ixo == root_sieve_size)
            ixo = 0;

    //assert(ixo == ix);
    pNV->root_idx = ixo;
    //pNH->current_value = val * (28ull/*root_sieve[0]*/ + 1);
    //pNV->root_idx = 1;
    //tpPrime cval = pNH->current_value;
    //tpRoot ix1 = pNV->root_idx;
    //while (cval < pNV->prime_value * pNV->prime_value)
    //{
    //    //compute new sieve value and advance Idx
    //    cval += pNV->prime_value * root_sieve[ix1++];
    //    //reset Idx if required
    //    if (ix1 == root_sieve_size)
    //        ix1 = 0;
    //}
    //assert(cval == pNV->prime_value * pNV->prime_value);
    //assert(ix == ix1);
}

void cHeapSieve::AddNewNode(tpPrime val)
{
    //advance and add a new row if necessary
    if (++idx_i == row_size)
    {
        AddNewRow();
        idx_i = 0;
    }

    InitializeNode(val_primes[idx_j] + idx_i, hep_primes[idx_j] + idx_i, val);
}

cHeapGenerator::cHeapGenerator(void) :
    writer("./PrimesOutput")
{
    InitializeRootSieve();
    InitializeHeapSieve();

    //get the future from the exit signal
    exitFuture = exitSignal.get_future();

    //manually add the root and its children to avoid 
    //complicated test for initial conditions in the loops
#ifdef USE_HEP_ROOT
    InitializeNode(&val_root, &hep_root, firstPrimes[root_N] + 0);
    for (auto i : Range<no_of_children>(1))
#else
    for (auto i : Range<manually_inserted_primes>(0))
#endif

    {
        AddNewNode(firstPrimes[root_N + i]);
        nr_hep_primes++;
    }
#ifdef USE_HEP_ROOT
    min_val_now = hep_root.current_value;
#else
    min_val_now = hep_primes[0][0].current_value;
#endif

    //account for first N+1 primes
    //N from root sieve + all manually added
    for (auto i : Range<root_N + manually_inserted_primes>())
    {
        no_gen_primes++;
        //checker.check_next_prime(firstPrimes[i]);
        assert(checker.check_next_prime(firstPrimes[i]));
    }
}

cHeapGenerator::~cHeapGenerator(void)
{
    DeinitializeRootSieve();
    DeinitializeHeapSieve();
}

void cHeapGenerator::PrintStatus(void)
{
    std::cout << std::setw(15) << candidate;
    std::cout << " [" << std::setw(2) << idx_j << "][";
    std::cout << std::setw(7) << idx_i << "]";
    std::cout << std::setw(16) << no_gen_primes << "\n";
}

void cHeapGenerator::Writer(std::future<void> futureObj)
{
    while (futureObj.wait_for(std::chrono::milliseconds(100)) == std::future_status::timeout)
    {
        WritePrimes(1);
    }
}

void cHeapGenerator::WritePrimes(int offset)
{
    static tpPrime cnt = 0;
    static tpRoot ix = 0;
    static tpRoot jx = 0;

    //capture the value at entry - it will be modified in Gen thread
    //but we only need to print some
    tpPrime n = nr_hep_primes - offset;

    while (cnt++ < n)
    {
#ifdef USE_STL_VECTOR
        resize_mutex.lock();
        //the STL vector is not thread safe
        //durign resize it will blow up
#endif
        tpPrime Nr{};// = vec_primes[jx][ix].prime_value; TO BE MODIFIED
#ifdef USE_STL_VECTOR
        resize_mutex.unlock();
#endif
        writer.WritePrime(Nr);
        //advance row if necessary
        if (++ix == row_size)
        {
            ix = 0;
            jx++;
        }
    }
}
