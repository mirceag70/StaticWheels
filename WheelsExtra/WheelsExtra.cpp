

#include <iostream>

void testHeapGenerator(const uint64_t LIMIT);
void test357(const uint64_t LIMIT);
void testCachedSieve(const uint64_t LIMIT, bool dogen);
void testCachedSieve2T(const uint64_t LIMIT);
void testAdvancedCachedSieve(const uint64_t LIMIT, bool dogen);
void testAdvancedCachedSieve2T(const uint64_t LIMIT);

constexpr uint64_t LIMIT = 1'000'000'000;
constexpr uint64_t LIMIT2 = 1ull<<38;

int main()
{
    std::locale mylocale("");   // get global locale 
    std::cout.imbue(mylocale);  // imbue global locale for thousands delimiter

    //testHeapGenerator(LIMIT);
    
    test357(LIMIT);
    testCachedSieve(LIMIT, true);
    testCachedSieve2T(LIMIT);
    testAdvancedCachedSieve(LIMIT, true);
    testAdvancedCachedSieve2T(LIMIT);

    testCachedSieve(LIMIT, false);
    testAdvancedCachedSieve(LIMIT2, false);
}

