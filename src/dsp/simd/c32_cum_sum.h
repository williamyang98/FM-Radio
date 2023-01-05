#pragma once

#include <immintrin.h>
#include <stdint.h>
#include <complex>
#include "simd_config.h"
#include "data_packing.h"

// Accumulate sum packed complex float
// NOTE: For performance we allow for modification of the input vector

#if defined(_DSP_SSSE3)
static inline 
std::complex<float> c32_cum_sum_ssse3(cpx128_t& a0) {
    // [c0 c1]
    cpx128_t a1;
    // [0 c0]
    a1.i = _mm_srli_si128(a0.i, 1*sizeof(std::complex<float>));
    // [0 c0+c1]
    a0.ps = _mm_add_ps(a0.ps, a1.ps);
    return a0.c32[0];
}
#endif

#if defined(_DSP_AVX2)
static inline 
std::complex<float> c32_cum_sum_avx2(cpx256_t& v_sum) {
    cpx128_t a0, a1;
    // [c0 c1]
    a0.ps = _mm_add_ps(v_sum.m128[0], v_sum.m128[1]);
    // [0 c0]
    a1.i = _mm_srli_si128(a0.i, 1*sizeof(std::complex<float>));
    // [0 c0+c1]
    a0.ps = _mm_add_ps(a0.ps, a1.ps);
    return a0.c32[0];
}
#endif