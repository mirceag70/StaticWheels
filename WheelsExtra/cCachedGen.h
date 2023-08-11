#pragma once
#include "generators.h"

class cCachedGen
{
private:
    tpSieveElement* sieve = nullptr;
    tpPrime sieve_size = 0;
    tpPrime limit;

    //tpPrime n_sum = 0;

    cChecker checker;
    cWriter writer;

    tpPrime CountSieve1(void);
    tpPrime CountSieve2(void);

    template<typename... Args>
    inline unsigned UnpackSlot(const tpSieveElement slot, tpPrime nk, const Args&&... rargs);

public:
    enum class UpperLimit : unsigned { OneBil = 3401, TenBil = 9592, HundredBil = 32000, ThousandBil = 78498 };

    tpPrime SoE_Rooted_30(void);
    tpPrime SoE_Rooted_Simple(void);

    cCachedGen(tpPrime lmt) : limit(lmt),
        writer("E:/list/out/PrimesOutput") {}
};

