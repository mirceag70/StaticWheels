#include "Helper.h"
//#define TESTWHEELSINITS
#include "wheels.h"

void test_initializations(void)
{
#ifdef TESTWHEELSINITS
    static_assert( 0 == wheels::ConstRound( 0));
    static_assert( 5 == wheels::ConstRound( 5.35));
    static_assert(-5 == wheels::ConstRound(-5.35));
    static_assert( 6 == wheels::ConstRound( 5.75));
    static_assert(-6 == wheels::ConstRound(-5.65));
    static_assert(wheels::ConstSqrt(25) == 5);
    static_assert(wheels::ConstSqrt(27345 * 27345) == 27345);
    static_assert(wheels::ConstSqrt(1'000'000'000'000) == 1'000'000);

    // test first primes
    for (auto i : wheels::range(wheels::_N))
        AddPrime(wheels::first97primes[i]);
    for (auto i : wheels::range(wheels::_N)) std::cout << wheels::first97primes[i] << "|";
    std::cout << "\n First primes OK.\n";

    wheels::RootSieve<3> rs; static_assert(rs.wheel_values == 8);
    for (auto i : wheels::range(rs.wheel_values)) std::cout << +rs.wheel[i] << "|";
    std::cout << " a total of " << rs.wheel_values << " values for wheel " << rs.root_N << ".\n Wheel OK.\n";

    wheels::RootSieve<5> rs5; static_assert(rs5.wheel_span == 2310);

    cTimer tmr; tmr.Start();
    //wheels::RootSieve<6> rs6; static_assert(rs6.wheel_span == 30030);
    tmr.Stop(true, "factors for 6"); nln();

    wheels::RootSieve<4> rs4;
    static_assert(rs4.root_N == 4);
    static_assert(rs4.wheel[4] == 19);
    for (auto i : wheels::range(rs4.wheel_values)) std::cout << +rs4.wheel[i] << "|";
    std::cout << " a total of " << rs4.wheel_values << " values for wheel " << rs4.root_N << ".\n Wheel OK.\n";

    for (auto i : wheels::range(rs.wheel_values))
    {
        for (auto j : wheels::range(rs.wheel_values))
            std::cout << +rs.factor[i][j] << "|";
        nln();
    }
    for (auto i : wheels::range(rs4.wheel_values))
    {
        for (auto j : wheels::range(rs4.wheel_values))
            std::cout << +rs4.factor[i][j] << "|";
        nln();
    }
    //for (auto i : wheels::range(rs5.wheel_values))
    //{
    //    for (auto j : wheels::range(rs5.wheel_values))
    //        std::cout << rs5.factor[i][j] << "|";
    //    nln();
    //}
#endif // TESTWHEELSINITS
}
