#define ITERATIONS	1

#include "Helper.h"

void InitializeRootSieve1(const unsigned char N);
void InitializeRootSieve2(const unsigned char N);
void InitializeRootSieve3(const unsigned char N);

uint64_t SoE_6k(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);
uint64_t SoE_rooted(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);
uint64_t SoE_rooted1(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);

uint64_t SoE_w8pat(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);

extern const unsigned rootN; 
const unsigned rootN = 4u;
constexpr uint64_t LIMIT = 1000000000;

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

int main()
{
    std::locale mylocale("");   // get global locale 
    std::cout.imbue(mylocale);  // imbue global locale for thousands delimiter

    //InitializeRootSieve1(rootN);
    //InitializeRootSieve2(rootN);
    //InitializeRootSieve3(rootN);

    //Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
    //    (LIMIT, "SoE 6k 1bit", &SoE_6k, false);
    //Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
    //    (LIMIT, "SoE rooted", &SoE_rooted, false);
    Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
        (LIMIT, "SoE rooted+", &SoE_rooted1, false);
    Try_Sieve<uint8_t, int, int, 8*(LIMIT / (8*30) + 1)>
        (LIMIT, "SoE pattern 8", &SoE_w8pat, false);

    return 0;
}
