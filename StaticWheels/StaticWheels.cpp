#ifdef DEBUG
#define ITERATIONS	1
#else
#define ITERATIONS	1
#endif // DEBUG


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
uint64_t SoE_bitpat_free0(const uint64_t Nmax, void*, void*, void*);
uint64_t SoE_bitpat_free1(const uint64_t Nmax, void*, void*, void*);
uint64_t SoE_pat_w7(const uint64_t Nmax, void*, void*, void*);

void InitBigWheel(void);

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
void test_wheels(void);
void testbigdiv(void);
void testalloc(void);

constexpr auto SZ = 1000u;
std::array<std::array<int, SZ>, SZ> mx;
void testtranspose(void)
{
    for (int i = 0; i < SZ; i++)
        for (int j = 0; j < SZ; j++)
            mx[i][j] = j;

    // cache friendly in place transposition
    constexpr unsigned block_size = 10; assert(SZ % block_size == 0);
    std::for_each(std::execution::par_unseq,
        range(SZ / block_size).begin(), range(SZ / block_size).end(), [&](const unsigned iib)
        {   //cache friendly transposition
            const unsigned ii = iib * block_size;
            //transpose the diagonal block
            for (unsigned i = ii; i < ii + block_size; i++)
                for (unsigned j = i + 1; j < ii + block_size; j++)
                    std::swap(mx[i][j], mx[j][i]);
            //transpose the rest of blocks
            for (unsigned jj = ii + block_size; jj < SZ; jj += block_size)
                for (unsigned i = ii; i < ii + block_size; i++)
                    for (unsigned j = jj; j < jj + block_size; j++)
                        std::swap(mx[i][j], mx[j][i]);
        });

    for (int i = 0; i < SZ; i++)
        for (int j = 0; j < SZ; j++)
            assert(mx[i][j] == i);
}

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


unsigned root_N = 4u;
constexpr uint64_t LIMIT = 1'000'000'000'000;
//constexpr uint64_t LIMIT = 0xffff'ffff;
//constexpr uint64_t LIMIT = (1ull << 36) - 0;
//constexpr uint64_t LIMIT = 100'00'000'00;

int main()
{
    std::locale mylocale("");   // get global locale 
    std::cout.imbue(mylocale);  // imbue global locale for thousands delimiter

    //testbigdiv(); return 0;
    //testalloc(); return 0;
    //testtranspose(); return 0;

    //testCSA(); return 0;
    //test_initializations(); return 0;
    //testforsw(4); return 0;
    //test_wheels(); return 0;

    //InitializeRootSieve1(root_N);         // brute force
    //InitializeRootSieve3(root_N);         // SoP
    //InitializeRootSieve2(root_N, true);   // 1bit 6k
    //*return 0;

    //Basic wheel based sieving
    //Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
    //    (LIMIT, "SoE 6k 1bit", &SoE_6k, false);
    //*Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
    //*    (LIMIT, "SoE rooted", &SoE_rooted, false);
    //Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
    //    (LIMIT, "SoE rooted+", &SoE_rooted1, false);
    
    //*Try_Sieve<uint8_t, int, int, 8*(LIMIT / (8*30) + 1)>
    //*    (LIMIT, "SoE pattern 8", &SoE_w8pat, false);
    //*Try_Sieve<uint8_t, int, int, 8 * (LIMIT / (8 * 30) + 1)>
    //*    (LIMIT, "SoE pattern 8 - incremental", &SoE_w8pat_inc, false);
    //*Try_Sieve<uint8_t, int, int, 8 * (LIMIT / (8 * 30) + 1)>
    //*    (LIMIT, "SoE pattern 8 - inc.vector", &SoE_w8pat_inc1, false);
    //W3 8T
    //Try_Sieve<uint8_t, int, int, 8 * (LIMIT / (8 * 30) + 1)>
    //    (LIMIT, "SoE pattern 8 - W3 8T", &SoE_w8pat_inc2, false);

    //Try_Sieve<int, int, int, 0>
    //    (LIMIT, "SoE pattern - free", &SoE_pat_free, false);   
    //*Try_Sieve<int, int, int, 0>
    //*    (LIMIT, "SoE pattern - bit", &SoE_bitpat_free0, false);
    //*Try_Sieve<int, int, int, 0>
    //*    (LIMIT, "SoE pattern - bit + list", &SoE_bitpat_free1, false);

    InitBigWheel();
    Try_Sieve<int, int, int, 0>
        (LIMIT, "SoE pattern - w7", &SoE_pat_w7, false);

    nln(); return 0;
}
