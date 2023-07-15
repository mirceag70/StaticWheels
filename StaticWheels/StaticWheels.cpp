#define ITERATIONS	1

#include "Helper.h"

void InitializeRootSieve1(const unsigned char N);
void InitializeRootSieve2(const unsigned char N, bool vb = false);
void InitializeRootSieve3(const unsigned char N);

uint64_t SoE_6k(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);
uint64_t SoE_rooted(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);
uint64_t SoE_rooted1(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);

uint64_t SoE_w8pat(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);
uint64_t SoE_w8pat_inc(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);
uint64_t SoE_w8pat_inc1(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);
uint64_t SoE_w8pat_inc2(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);

uint64_t SoE_pat_free(const uint64_t Nmax, void*, void*, void*);

unsigned root_N = 7u;
//constexpr uint64_t LIMIT = 1'000'000'000'000;
//constexpr uint64_t LIMIT = 0xffff'ffff;
//constexpr uint64_t LIMIT = (1ull << 32) - 0;
constexpr uint64_t LIMIT = 100'000'000'000'000;

//uint64_t idx2no30_armtc(uint64_t idx)
//{
//    uint64_t k = idx / 8;   uint64_t i = idx % 8;
//    return 4 * (i + (i / 6)) + 7 - 2 * ((i / 2) + (i / 7)) + 30 * k;
//}
//uint64_t idx2no30_logic(uint64_t idx)
//{
//    uint64_t k = idx / 8;   uint64_t i = idx % 8;
//    switch (i)
//    {
//    case 0: i = 7; break; case 1: i = 11; break;
//    case 2: i = 13; break; case 3: i = 17; break;
//    case 4: i = 19; break; case 5: i = 23; break;
//    case 6: i = 29; break; case 7: i = 31; break;
//    default: assert(false);
//    }
//    return i + 30 * k;
//}
//uint64_t no2idx30_armtc(uint64_t idx)
//{
//    uint64_t k = idx / 8;   uint64_t i = idx % 8;
//    return 4 * (i + (i / 6)) + 7 - 2 * ((i / 2) + (i / 7)) + 30 * k;
//}
//uint64_t no2idx30_logic(uint64_t idx)
//{
//    uint64_t k = idx / 8;   uint64_t i = idx % 8;
//    switch (i)
//    {
//    case 0: i = 7; break; case 1: i = 11; break;
//    case 2: i = 13; break; case 3: i = 17; break;
//    case 4: i = 19; break; case 5: i = 23; break;
//    case 6: i = 29; break; case 7: i = 31; break;
//    default: assert(false);
//    }
//    return i + 30 * k;
//}
//void testfn(void)
//{
//    cTimer tmr; tmr.Start();
//
//    uint64_t c = 0;
//
//    for (uint64_t i = 0; i < 10'000'000'000; i++)
//    {
//        uint64_t n = idx2no30_armtc(i);
//        //uint64_t n1 = idx2no30_logic(i);
//        //assert(n == n1);
//        if ((n & 1) == 1)
//        {
//            //std::cout << "!!";
//            c++;
//        }
//        //switch (k)
//        //{
//        //default: assert(false);
//        //case 1: case 7: case 11: case 13:
//        //case 17: case 19: case 23: case 29:
//        //    ;
//        //}
//    }
//    std::cout << c;
//    tmr.Stop(true);
//}

void testCSA(void);
void test_initializations(void);
void testforsw(unsigned start)
{
    constexpr unsigned sz = 100;
    unsigned arr[sz]{};
    unsigned i = start;

    switch (start)
    {
        default:    NEVERHERE;

        for (;;)
        {
            case 1: arr[i++] = 55; if (i == sz) break;
            case 2: arr[i++] = 55; if (i == sz) break;
            case 3: arr[i++] = 55; if (i == sz) break;
            case 4: arr[i++] = 55; if (i == sz) break;
            case 5: arr[i++] = 55; if (i == sz) break;
        }
    }
}

int main()
{
    std::locale mylocale("");   // get global locale 
    std::cout.imbue(mylocale);  // imbue global locale for thousands delimiter

    //testCSA(); return 0;
    //test_initializations(); return 0;
    //testforsw(4); return 0;

    //InitializeRootSieve1(root_N);
    //InitializeRootSieve2(root_N); return 0;
    //InitializeRootSieve3(root_N);

    //Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
    //    (LIMIT, "SoE 6k 1bit", &SoE_6k, false);
    //Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
    //    (LIMIT, "SoE rooted", &SoE_rooted, false);
    //Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
    //    (LIMIT, "SoE rooted+", &SoE_rooted1, false);
    //Try_Sieve<uint8_t, int, int, 8*(LIMIT / (8*30) + 1)>
    //    (LIMIT, "SoE pattern 8", &SoE_w8pat, false);
    //Try_Sieve<uint8_t, int, int, 8 * (LIMIT / (8 * 30) + 1)>
    //    (LIMIT, "SoE pattern 8 - incremental", &SoE_w8pat_inc, false);
    //Try_Sieve<uint8_t, int, int, 8 * (LIMIT / (8 * 30) + 1)>
    //    (LIMIT, "SoE pattern 8 - inc.vector", &SoE_w8pat_inc1, false);
    //Try_Sieve<uint8_t, int, int, 8 * (LIMIT / (8 * 30) + 1)>
    //    (LIMIT, "SoE pattern 8 - inc.vector.full", &SoE_w8pat_inc2, false);
    Try_Sieve<int, int, int, 0>
        (LIMIT, "SoE pattern - free", &SoE_pat_free, false);
    
    return 0;
}
