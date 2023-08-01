#pragma once

#include <immintrin.h>
#include <stdint.h>
#include <stdalign.h>
#include "../simd_flags.h"

// Accumulate sum packed float
#if defined(__SSE2__)
static inline 
float f32_cum_sum_sse2(__m128 x) {
    // [f0 f1 f2 f3]
    __m128i a0 = _mm_castps_si128(x);
    // [f2 f3 0 0]
    __m128i a1 = _mm_srli_si128(a0, 2*sizeof(float));
    // [f0+f2 f1+f3 0 0]
    __m128 b0 = _mm_add_ps(_mm_castsi128_ps(a0), _mm_castsi128_ps(a1));
    __m128i a2 = _mm_castps_si128(b0);
    // [f1+f3 0 0 0]
    __m128i a3 = _mm_srli_si128(a2, 1*sizeof(float));
    // [f0+f1+f2+f3 0 0 0]
    __m128 b1 = _mm_add_ps(_mm_castsi128_ps(a2), _mm_castsi128_ps(a3));

    alignas(16) float y[4];
    _mm_store_ps(y, b1);
    return y[0];
}
#endif

#if defined(__AVX__)
static inline 
float f32_cum_sum_avx(__m256 x) {
    // [f0+f4 f1+f5 f2+f6 f3+f7]
    __m128 a0 = _mm_add_ps(
        _mm256_extractf128_ps(x, 0),
        _mm256_extractf128_ps(x, 1)
    );
    return f32_cum_sum_sse2(a0);
}
#endif