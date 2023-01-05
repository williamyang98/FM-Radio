#pragma once

#include <immintrin.h>
#include <stdint.h>
#include "simd_config.h"
#include "data_packing.h"

// Accumulate sum packed float
// NOTE: For performance we allow for modification of the input vector

#if defined(_DSP_SSSE3)
static inline 
float f32_cum_sum_ssse3(cpx128_t& a0) {
    // [f0 f1 f2 f3]
    cpx128_t a1;
    // [0 0 f0 f1]
    a1.i = _mm_srli_si128(a0.i, 2*sizeof(float));
    // [0 0 f0+f2 f1+f3]
    a0.ps = _mm_add_ps(a0.ps, a1.ps);
    // [0 0 0 f0+f2]
    a1.i = _mm_srli_si128(a0.i, 1*sizeof(float));
    // [0 0 0 f0+f1+f2+f3]
    a0.ps = _mm_add_ps(a0.ps, a1.ps);
    return a0.f32[0];
}
#endif

#if defined(_DSP_AVX2)
static inline 
float f32_cum_sum_avx2(cpx256_t& v_sum) {
    cpx128_t a0, a1;
    // [f0 f1 f2 f3]
    a0.ps = _mm_add_ps(v_sum.m128[0], v_sum.m128[1]);
    // [0 0 f0 f1]
    a1.i = _mm_srli_si128(a0.i, 2*sizeof(float));
    // [0 0 f0+f2 f1+f3]
    a0.ps = _mm_add_ps(a0.ps, a1.ps);
    // [0 0 0 f0+f2]
    a1.i = _mm_srli_si128(a0.i, 1*sizeof(float));
    // [0 0 0 f0+f1+f2+f3]
    a0.ps = _mm_add_ps(a0.ps, a1.ps);
    return a0.f32[0];
}
#endif