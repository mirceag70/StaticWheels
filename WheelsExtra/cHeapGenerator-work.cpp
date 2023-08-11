#include "cHeapGenerator.h"

constexpr auto WAIT_TIMEOUT = std::chrono::seconds(10);

tpPrime cHeapGenerator::GeneratePrimes(tpPrime limit)
{
    //because I/O ops are time consuming
    //the primes are to be generated in a parallel thread
    // for which we need to create a signal to stop it:
    std::promise<void> exitSignalW;
    std::future<void> futureObjW = exitSignalW.get_future();

    cTimer tmr;
    tmr.Start();

    std::thread tg(&cHeapGenerator::Generator, this, limit);
    //std::thread tw(&cHeapGenerator::Writer, this, std::move(futureObjW));

    while (exitFuture.wait_for(WAIT_TIMEOUT) == std::future_status::timeout)
    {
        PrintStatus();
    }

    tg.join();
    exitSignalW.set_value();
    //tw.join();

    //write down the rest of the primes, if any
    //WritePrimes();

    tmr.Stop(true); nln();

    PrintStatus();

    return no_gen_primes;
}

unsigned cnt0[] = { 0, 0, 0, 0 };

void cHeapGenerator::Generator(tpPrime lmt)
{
    limit = lmt;
    sqrtlimit = (tpPrime)(sqrt(limit)+1);

    cTimer tmr;

    tmr.Start();

    tpRoot i = 0;
    //some primes already placed in sieve
    //so we must skip them 
    for (; i  < manually_inserted_primes+1; i++)
        candidate += root_sieve[i];

    do
    {
        Process();
        candidate += root_sieve[i++];
        if (i == root_sieve_size)
            i = 0;
    } 
    while (candidate < limit);

    tmr.Stop(true); nln();
    exitSignal.set_value();

    //std::cout << "Rsd:" << cnt_rsd << " - Eg:" << cnt_eq << " - Lw:"; //nln();
    //std::cout << "Adv:" << cnt_adv << " - Swc:" << cnt_swtchright; nln();
    //std::cout << "0:" << cnt0[0] << " 1:" << cnt0[1] << " 2:" << cnt0[2] << " 3:" << cnt0[3];
}

void cHeapGenerator::Process(void)
{
    if (candidate == min_val_now)
    {
        AdvanceSieve();
        //cnt_eq++;
    }
    else
    {
        assert(candidate < min_val_now);

        InsertNewPrime();
    }
}

void cHeapGenerator::InsertNewPrime(void)
{
    //check value
    //if(!checker.check_next_prime(candidate))
    //    __debugbreak();;
    assert(checker.check_next_prime(candidate));
    no_gen_primes++;
    if (candidate <= sqrtlimit)
    {
        AddNewNode(candidate);
        nr_hep_primes++;
    }
}

//__declspec(noinline)
void cHeapGenerator::AdvanceSieve(void)
{
#ifdef USE_HEP_ROOT
    cPrimeHep* pF = &hep_root;
#else
    cPrime* pF = hep_primes[0] + 0;
#endif
    
    // advance while root == min val
    do
    {
        //compute new sieve value and advance Idx
        pF->current_value += pF->prime->prime_value * root_sieve[pF->prime->root_idx++];
        //reset Idx if required
        if (pF->prime->root_idx == root_sieve_size)
            pF->prime->root_idx = 0;

        Heapify4(*pF);

        //cnt_adv++;
    } 
    while (pF->current_value == min_val_now);

    min_val_now = pF->current_value;
}

/*
void cHeapGenerator::Heapify(void)
{
    tpPrime valnow = hep_root.current_value;

    cPrime* pChild = hep_primes[0] + 0;
    cPrime* pChild_R = pChild + 1;

    tpPrime idxi_current = 0;
    //first we deal with the root and its children
    if (pChild->current_value <= pChild_R->current_value)
        idxi_current = 0;
    else
    {
        idxi_current = 1;
        pChild = pChild_R;
    }

    if (pChild->current_value >= valnow)
        //nothing to do
        return;

    tpPrime idxj_current = 0;

    cPrime OldRoot = hep_root;
    hep_root = *pChild;

    tpPrime last_idx_j = idx_j;
    //go back one place if not R child
    tpPrime last_idx_i = (idx_i % 2 == 1) ? idx_i : idx_i - 1;

    // work down from the next first child of idx_smlst
    tpPrime idxi_child = (idxi_current + 1) * no_of_children;
    tpPrime idxj_child = 0;

    //move down nodes until a larger node
    do
    {
        // we must choose the lowest child to bubble it up
        pChild = hep_primes[idxj_child]+idxi_child;
        pChild_R = pChild+1;

        //fixate on the smallest child
        if (pChild->current_value > pChild_R->current_value)
        {
            idxi_child++;
            pChild = pChild_R;

        }
        //if chosen child is lower, raise it
        if (pChild->current_value < valnow)
        {
            hep_primes[idxj_current][idxi_current] = *pChild;
            cnt_rsd++;
        }
        else
        {
            //we have reached a higher or equal value, so stop here
            break;
        }

        //advance the chosen child as current parrent 
        idxi_current = idxi_child;
        idxj_current = idxj_child;

        // new idx: 2*(i+1) ; IS NOT 2i+1 BECAUSE ROOT IS NOT IN HEAP!
        idxi_child++;
        idxi_child <<= 1;
        idxj_child <<= 1;
        while (idxi_child >= row_size)
        {
            idxi_child -= row_size;
            idxj_child++;
        }
        assert(idxi_child < row_size);

        //as long as we have a child, continue
        if (idxj_child < last_idx_j)
            continue;
        else
            if (idxj_child > last_idx_j)
                break;
            else
                if (idxi_child > last_idx_i)
                    break;
    } while (true);

    //and finally place old root it where it belongs
    hep_primes[idxj_current][idxi_current] = OldRoot;
}
*/


//extern "C" tpPrime HeapifyASM(int, int, cPrime***, cPrime*);

//__declspec(noinline)
void cHeapGenerator::Heapify2(cPrimeHep OldRoot)
{
    // work down from the root
#ifdef USE_HEP_ROOT
    cPrimeHep* pCurrent = &hep_root;
    tpPrime idxi_child = 0;
#else
    cPrime* pCurrent = hep_primes[0] + 0;
    tpPrime idxi_child = 1;
#endif
    
    //pRoot is the root
    tpPrime valnow = OldRoot.current_value;

    tpPrime idxj_child = 0;

    //if last leaf has no right sibling, ignore it and consider sieves stops at the previous leaf!
    //(anyway, it is very, very big relative to the others, we'll never have to swap it)
    //thus we avoid one complex test
    tpPrime child_ix = idxi_child;
    tpPrime last_ix = nr_hep_primes-1;
    if (last_ix % 2 == 1)
        last_ix--;

    //move down nodes until a larger node
    do
    {
        cPrimeHep* pChild = hep_primes[idxj_child] + idxi_child;

        // we must choose the lowest child to bubble it up
#ifdef USE_HEP_ROOT
        //idxi_childR = idxi_child + 1;
        //idxj_childR = idxj_child;
        //cPrime* pChild_R = hep_primes[idxj_childR] + idxi_childR;
        cPrimeHep* pChild_R = pChild+1;
        //we must choose R if lower than L
        if (pChild_R->current_value < pChild->current_value)
        {
            idxi_child++;// = idxi_childR;
            //idxj_child = idxj_childR;
            pChild = pChild_R;
            child_ix++;
        }
#else
        tpPrime idxi_childR, idxj_childR;
        //if we are at the end of a row, go to next row
        if ((idxi_childR = idxi_child + 1) < row_size)
        {
            idxj_childR = idxj_child;
        }
        else
        {
            idxi_childR = 0;
            idxj_childR = idxj_child + 1;
        }
        cPrime* pChild_R = hep_primes[idxj_childR] + idxi_childR;

        //we must choose R if lower than L
        if (pChild_R->current_value < pChild->current_value)
        {
            idxi_child = idxi_childR;
            idxj_child = idxj_childR;
            pChild = pChild_R;
            child_ix++;
        }
#endif

        //if chosen child is lower, raise it
        if (pChild->current_value < valnow)
        {
            *pCurrent = *pChild;
            //cnt_rsd++;
        }
        else
        {
            //we have reached a higher or equal value, so stop here
            break;
        }

        //advance the chosen child as current parrent 
        pCurrent = pChild;

#ifdef USE_HEP_ROOT
        // new idx: 2*(i+1) ; NOT 2i+1 BECAUSE ROOT IS NOT IN HEAP!
        idxi_child++;
        idxi_child <<= 1;
        idxj_child <<= 1;
        if(idxi_child >= row_size)
        {
            idxi_child -= row_size;
            idxj_child++;
            if (idxi_child >= row_size)
            {
                idxi_child -= row_size;
                idxj_child++;
            }
        }
        child_ix++;
        child_ix *= 2;
        assert(child_ix == idxj_child * row_size + idxi_child);
#else
        // new idx: 2*i + 1 
        idxi_child <<= 1;
        idxj_child <<= 1;
        if ((++idxi_child) >= row_size)
        {
            idxi_child -= row_size;
            idxj_child++;
        }
        child_ix *= 2;
        child_ix++;
        assert(child_ix == idxj_child * row_size + idxi_child);
#endif

        //as long as we have a child, continue
        if(child_ix > last_ix)
            break;

    } while (true);

    //and finally place it where it belongs
    *pCurrent = OldRoot;
}

//__declspec(noinline)
void cHeapGenerator::AdvanceIdxs4(tpPrime& idxi_child, tpPrime& idxj_child, const tpPrime idx_smlst)
{
    // new idx: 4*(i+1) ; NOT 4i+1 BECAUSE ROOT IS NOT IN HEAP!
    idxi_child += idx_smlst + 1;
    idxi_child <<= 2;
    idxj_child <<= 2;
    while (idxi_child >= row_size)
    {
        idxi_child -= row_size;
        idxj_child++;
    }
}

//__declspec(noinline)
unsigned SmlstOf4Idx(tpPrime a, tpPrime b, tpPrime c, tpPrime d)
{
    if (a <= b)     //choose between a, c, d
        if (a <= c) //choose between a, d
            if (a <= d) return 0;
            else        return 3;
        else        //choose between c, d
            if (c <= d) return 2;
            else        return 3;
    else            //choose between b, c, d
        if (b <= c) //choose between b, d
            if (b <= d) return 1;
            else        return 3;
        else        //choose between c, d
            if (c <= d) return 2;
            else        return 3;
}

unsigned SmlstOf4Idx(cPrimeHep* pP)
{
    unsigned ix1 = 1 - (pP[0].current_value <= pP[1].current_value);
    unsigned ix2 = 3 - (pP[2].current_value <= pP[3].current_value);
    return (pP[ix1].current_value <= pP[ix2].current_value) ? ix1 : ix2;
    // 15% worse:
    //unsigned a = pP[ix1].current_value <= pP[ix2].current_value;
    //return a * ix1 + (1 - a) * ix2;
}

//__declspec(noinline)
void cHeapGenerator::Heapify4(cPrimeHep OldRoot)
{
    // work down from the root
    //tpPrime a = hep_primes[0][0].current_value;
    cPrimeHep* pCurrent = &hep_root;
    tpPrime idxi_child = 0;
    tpPrime idxj_child = 0;

    tpPrime valnow = OldRoot.current_value;

    //break test indexes 
    tpPrime child_ix = idxi_child;
    tpPrime last_ix = nr_hep_primes;
    //if last set of children is not complete, advance to the end of its space
    while ((last_ix % no_of_children) != 0) last_ix++;
    last_ix--;

    //move down nodes until a larger node
    //do
    for(;;)
    {
        cPrimeHep* pChild = hep_primes[idxj_child] + idxi_child;
        //tpPrime idx_smlst = SmlstOf4Idx( //a,
        //                        (pChild + 0)->current_value,
        //                        (pChild + 1)->current_value,
        //                        (pChild + 2)->current_value,
        //                        (pChild + 3)->current_value);
        tpPrime idx_smlst = SmlstOf4Idx(pChild);
        //cnt0[idx_smlst]++;
        //unsigned idx_smlst = SmlstOf4Idx(pChild);
        //AdvanceIdxs4(idxi_child, idxj_child, idx_smlst);
        //a = hep_primes[idxj_child][idxi_child].current_value;

        // we must choose the lowest child to bubble it up
        pChild += idx_smlst;

        //if chosen child is lower, raise it
        if (pChild->current_value < valnow)
        //if (pChild->current_value < OldRoot.current_value)
        {
            *pCurrent = *pChild;
            //cnt_rsd++;
        }
        else
        {
            //we have reached a higher or equal value, so stop here
            break;
        }

        //advance the chosen child as current parrent 
        pCurrent = pChild;

        // new idx: 4*(i+1) ; NOT 4i+1 BECAUSE ROOT IS NOT IN HEAP!
        child_ix += idx_smlst + 1;
        idxi_child += idx_smlst + 1;
        idxj_child <<= 2;
        child_ix <<= 2;
        //assert(child_ix == idxj_child * row_size + idxi_child);
        idxi_child <<= 2;
        while (idxi_child >= row_size)
        {
            idxi_child -= row_size;
            idxj_child++;
        }

        //as long as we have a child, continue
        if (child_ix > last_ix)
            break;

    } //while (true);

    //and finally place it where it belongs
    *pCurrent = OldRoot;
}

void testHeapGenerator(const uint64_t LIMIT)
{
    cTimer tmr;
    tmr.Start();

    std::vector<decltype(tmr.LapTime())> times;

    for (auto i : Range<1>(1))
    {
        cHeapGenerator worker;
        //worker.GeneratePrimes(LIMIT); nln();
        worker.Generator(LIMIT); nln();
        std::cout << worker.no_gen_primes << " primes up to " << LIMIT; nln();
        times.push_back(tmr.LapTime(true));
        nln(true);
    }

    tmr.Stop();
    nln(true);
    std::cout << "Average compute time: " << Average(times);
    nln(true);
}

/*
//__declspec(noinline)
void cHeapGenerator::Heapify(cPrime OldRoot)
{
    //HeapifyASM(idx_j, idx_i, hep_primes, pRoot);
    //return;

    // work down from the root
    tpPrime idxi_current = 0;
    tpPrime idxj_current = 0;
    //pRoot is the root
    tpPrime valnow = OldRoot.current_value;

    //we know the first two children exists - we inserted three primes from the start
    tpPrime idxi_child = 1;
    tpPrime idxj_child = 0;

    tpPrime last_idx_j = idx_j;
    tpPrime last_idx_i = idx_i;
    //tpPrime last_idx_i = idx_i & ~1;
    //if last leaf has no right sibling, ignore it and consider sieves stops at the previous leaf!
    //(anyway, it is very, very big relative to the others, we'll never have to swap it)
    //thus we avoid one complex test
    if (last_idx_i % 2 == 1)
        last_idx_i--;
    //last_idx_i -= last_idx_i % 2;

    //move down nodes until a larger node
    do
    {
        tpPrime idxi_childR, idxj_childR;

        // we must choose the lowest child to bubble it up
        cPrime* pChild = hep_primes[idxj_child] + idxi_child;
        cPrime* pChild_R = pChild+1;

        //if we are at the end of a row, go to next row
        if ((idxi_childR = idxi_child + 1) < row_size)
        {
            idxj_childR = idxj_child;
        }
        else
        {
            idxi_childR = 0;
            idxj_childR = idxj_child + 1;
            //pChild_R = hep_primes[idxj_childR] + idxi_childR;
        }

        pChild_R = hep_primes[idxj_childR] + idxi_childR;

        //we must choose R if lower than L
        if (pChild_R->current_value < pChild->current_value)
        {
            idxi_child = idxi_childR;
            idxj_child = idxj_childR;
            pChild = pChild_R;

            //cnt_swtchright++;
        }
        //bool needSwitch = pChild_R->current_value < pChild->current_value;
        //idxi_child = needSwitch * idxi_childR + (true - needSwitch) * idxi_child;
        //idxj_child = needSwitch * idxj_childR + (true - needSwitch) * idxj_child;
        //pChild = (cPrime*)(needSwitch * (size_t)pChild_R + (true - needSwitch) * (size_t)pChild);

        //if chosen child is lower, raise it
        if (pChild->current_value < valnow)
        {
            hep_primes[idxj_current][idxi_current] = *pChild;

            //cnt_rsd++;
        }
        else
        {
            //we have reached a higher or equal value, so stop here
            break;
        }

        //advance the chosen child as current parrent
        idxi_current = idxi_child;
        idxj_current = idxj_child;

        // new idx: 2*i + 1
        idxi_child <<= 1;
        idxj_child <<= 1;
        if ((++idxi_child) >= row_size)
        {
            idxi_child -= row_size;
            idxj_child++;
        }

        //as long as we have a child, continue
        if (idxj_child < last_idx_j)
            continue;
        else
            if (idxj_child > last_idx_j)
                break;
            else
                if (idxi_child > last_idx_i)
                    break;
    } while (true);

    //and finally place it where it belongs
    hep_primes[idxj_current][idxi_current] = OldRoot;
}
*/