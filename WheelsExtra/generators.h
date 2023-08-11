#pragma once
#include "Helper.h"

constexpr auto CACHE_LINE_BYTES = 32;
#define CACHE_ALIGN __declspec(align(CACHE_LINE_BYTES))

typedef  uint8_t tpSieveElement;
static_assert(sizeof(tpSieveElement) == 1);

constexpr unsigned root_N = 3;
constexpr unsigned root_sieve_size = 8;
constexpr unsigned root_primorial = 2 * 3 * 5;
constexpr CACHE_ALIGN std::array<tpSieveElement, root_sieve_size> root_sieve_values
               { 1, 7, 11, 13, 17, 19, 23, 29};
constexpr CACHE_ALIGN std::array<tpSieveElement, root_sieve_size> root_sieve_steps
               { 6, 4, 2, 4, 2, 4, 6, 2 };

consteval uint8_t Msk(unsigned i)
{
    return (1 << (i));
};

constexpr CACHE_ALIGN std::array<tpSieveElement, root_primorial> root_sieve_masks
               { 0, Msk(0), 0, 0, 0,
                 0, 0, Msk(1), 0, 0,
                 0, Msk(2), 0, Msk(3), 0,
                 0, 0, Msk(4), 0, Msk(5),
                 0, 0, 0, Msk(6), 0,
                 0, 0, 0, 0, Msk(7) };
               //{200, Msk(0), 202, 203, 204,
               // 205, 206, Msk(1), 208, 209,
               // 210, Msk(2), 212, Msk(3), 214,
               // 215, 216, Msk(4), 218, Msk(5),
               // 220, 221, 222, Msk(6), 224,
               // 225, 226, 227, 228, Msk(7) };

constexpr CACHE_ALIGN std::array<tpSieveElement, root_primorial> bit_masks
        { Msk(0), Msk(1), Msk(2), Msk(3), 
          Msk(4), Msk(5), Msk(6), Msk(7) };

template<typename T1, typename T2>
constexpr bool GETbt(T1& vect, const T2 pos)
{
    T2 q = pos / 30;
    T2 r = pos % 30;
    tpSieveElement msk = root_sieve_masks[r];
    //assert(msk < 200);
    return (vect[q] & msk);
}

template<typename T1, typename T2>
constexpr void RESETbt(T1& vect, const T2 pos)
{
    T2 q = pos / 30;
    T2 r = pos % 30;
    tpSieveElement msk = root_sieve_masks[r];
    //assert(msk < 200);
    vect[q] &= ~msk;
}

constexpr void IncrementRootIdx(uint8_t& idx)
{
    idx == (root_sieve_size - 1) ? idx = 0 : idx++;
}