#pragma once

#include <immintrin.h>
#include <stdint.h>
#include <complex>
#include <stdalign.h>
#include "../simd_flags.h"

// Accumulate sum packed complex float
#if defined(__SSE2__)
static inline
std::complex<float> c32_cum_sum_sse2(__m128 x) {
    // [c0 c1]
    __m128i a0 = _mm_castps_si128(x);
    // [0 c0]
    __m128i a1 = _mm_srli_si128(a0, 1*sizeof(std::complex<float>));
    // [0 c0+c1]
    __m128 b0 = _mm_add_ps(
        _mm_castsi128_ps(a0),
        _mm_castsi128_ps(a1)
    );

    alignas(16) std::complex<float> y[2];
    _mm_store_ps(reinterpret_cast<float *>(y), b0);
    return y[0];
}
#endif

#if defined(__AVX__)
static inline
std::complex<float> c32_cum_sum_avx(__m256 x) {
    // [c0+c2 c1+c3]
    __m128 a0 = _mm_add_ps(
        _mm256_extractf128_ps(x, 0),
        _mm256_extractf128_ps(x, 1)
    );
    return c32_cum_sum_sse2(a0);
}
#endif