#define ITERATIONS	5

#include "Helper.h"

void InitializeRootSieve1(const unsigned char N);
void InitializeRootSieve2(const unsigned char N);
void InitializeRootSieve3(const unsigned char N);

uint64_t SoE_6k(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);
uint64_t SoE_rooted(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);
uint64_t SoE_rooted1(const uint64_t Nmax, uint8_t vPrimes[], void*, void*);

extern const unsigned rootN; 
const unsigned rootN = 9u;
constexpr uint64_t LIMIT = 1000000000;

int main()
{
    std::locale mylocale("");   // get global locale 
    std::cout.imbue(mylocale);  // imbue global locale for thousands delimiter

    //InitializeRootSieve1(rootN);
    //InitializeRootSieve2(rootN);
    //InitializeRootSieve3(rootN);

    Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
        (LIMIT, "SoE 6k 1bit", &SoE_6k, false);
    Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
        (LIMIT, "SoE rooted", &SoE_rooted, false);
    Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
        (LIMIT, "SoE rooted+", &SoE_rooted1, false);

    return 0;
}
