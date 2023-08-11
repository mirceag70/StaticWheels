#include "helper.h"

constexpr uint64_t szM = 450'000; bool vM[szM + 1];
uint64_t Sieve357(const uint64_t Nmax,
    uint64_t IP[], uint64_t IQ[], int32_t JQ[], bool initialize = false)
{
    uint64_t numPrimes = 0;     //for this call
    for (int i = 0; i < szM; i++) vM[i] = true;

    static uint64_t nstart, nend;
    static uint32_t nq, stepq;
    static int32_t iq;      //idx of the last root prime generated so far
    if (initialize)
    {
        vM[1] = false;
        JQ[0] = 3;  IQ[0] = 9;
        iq = 0; stepq = 2; nq = 5;
        nstart = 0; nend = szM;
        numPrimes = 1; AddPrime(2); // account for 2
    }

    uint64_t n, nsqrt = (uint64_t)ceil(sqrt(std::min(Nmax, nstart + szM)));
    //generate rest of root primes for this iteration
    for (; nq <= nsqrt; nq += stepq, stepq = 6 - stepq)
    {
        bool isprime = true;
        for (int i = 2; i <= iq; i++)
        {   //trial-division
            if ((nq / JQ[i]) * JQ[i] == nq)
            {
                isprime = false;
                break;
            }
        }
        if (isprime)
        {
            iq++;
            JQ[iq] = nq;
            IQ[iq] = ((uint64_t)nq) * nq;
        }
    }
    //strike out all composites
    for (int i = 0; i <= iq; i++)
    {
        uint32_t stp = 2 * JQ[i];
        for (n = IQ[i]; n < nend; n += stp)
        {
            assert(n - nstart < szM);
            vM[n - nstart] = false;
        }
        IQ[i] = n;
    }
    //move primes in IP
    for (int i = 1; i < szM; i += 2)
        if ((i + nstart) > Nmax)
            break;
        else
            if (vM[i])
                AddPrime(IP[numPrimes++] = nstart + i);
    //prepare for next call
    nstart = nend; nend += szM;

    return numPrimes;
}

void test357(const uint64_t LIMIT)
{
    std::cout << "\n - Sieve 357 - \n";
    uint64_t szall = PI_Nmax(LIMIT);
    uint64_t* vulPrimes = new uint64_t[szall + 1];
    uint64_t szroot = PI_Nmax((uint64_t)ceil(sqrt(std::max(LIMIT, szM))));
    //std::cout << sz;
    uint64_t* IQ = new uint64_t[szroot + 1];
    int32_t* JQ = new int32_t[szroot + 1];

    cTimer tmr;
    std::vector<decltype(tmr.LapTime())> times;

    tmr.Start();
    
    for (auto i : Range<1>())
    {
        nln();

        uint64_t numPrimes = Sieve357(LIMIT, vulPrimes, IQ, JQ, true);
        for (uint64_t start = szM; start < LIMIT; start += szM)
            numPrimes += Sieve357(LIMIT, vulPrimes, IQ, JQ);

        std::cout << numPrimes << " primes up to " << LIMIT;
        times.push_back(tmr.LapTime(true));
    }
    nln(true);
    tmr.Stop();
    std::cout << "Average compute time: " << Average(times);

    nln(true);
    times.clear();

    delete[] vulPrimes;
    delete[] IQ;
    delete[] JQ;
}