#pragma once
#include <assert.h>

// NOTE: Assumes arrays are aligned
// Multiply and accumulate vector of floats with another vector of floats

static inline
float f32_cum_mul_scalar(const float* x0, const float* x1, const int N) {
    float y = 0;
    for (int i = 0; i < N; i++) {
        y += x0[i] * x1[i];
    }
    return y;
}

// TODO: Modify code to support ARM platforms like Raspberry PI using NEON
#include <immintrin.h>
#include "simd_config.h"
#include "data_packing.h"

#if defined(_DSP_SSSE3)
static inline
float f32_cum_mul_ssse3(const float* x0, const float* x1, const int N)
{
    float y = 0;

    // 128bits = 16bytes = 4*4bytes
    const int K = 4;
    const int M = N/K;

    for (int i = 0; i < M; i++) {
        __m128 a0 = _mm_load_ps(&x0[i*K]);
        __m128 a1 = _mm_load_ps(&x1[i*K]);
        cpx128_t b0;
        b0.ps = _mm_mul_ps(a0, a1);

        // Perform vectorised cumulative sum
        // Shift half of vector and add. Repeat until we get the final sum
        cpx128_t b1;

        b1.i = _mm_srli_si128(b0.i, 2*sizeof(float));
        b0.ps = _mm_add_ps(b0.ps, b1.ps);
        b1.i = _mm_srli_si128(b0.i, 1*sizeof(float));
        b0.ps = _mm_add_ps(b0.ps, b1.ps);

        y += b0.f32[0];
    }

    const int N_vector = M*K;
    const int N_remain = N-N_vector;
    y += f32_cum_mul_scalar(&x0[N_vector], &x1[N_vector], N_remain);

    return y;
}
#endif

#if defined(_DSP_AVX2)
static inline
float f32_cum_mul_avx2(const float* x0, const float* x1, const int N)
{
    float y = 0;

    // 256bits = 32bytes = 8*4bytes
    const int K = 8;
    const int M = N/K;

    for (int i = 0; i < M; i++) {
        __m256 a0 = _mm256_load_ps(&x0[i*K]);
        __m256 a1 = _mm256_load_ps(&x1[i*K]);
        cpx256_t b0;
        b0.ps = _mm256_mul_ps(a0, a1);

        // Perform vectorised cumulative sum
        // Shift half of vector and add. Repeat until we get the final sum
        cpx128_t c0;
        cpx128_t c1;
        // merge 128bit lanes
        c0.ps = _mm_add_ps(b0.m128[0], b0.m128[1]);

        c1.i = _mm_srli_si128(c0.i, 2*sizeof(float));
        c0.ps = _mm_add_ps(c0.ps, c1.ps);
        c1.i = _mm_srli_si128(c0.i, 1*sizeof(float));
        c0.ps = _mm_add_ps(c0.ps, c1.ps);

        y += c0.f32[0];
    }

    const int N_vector = M*K;
    const int N_remain = N-N_vector;
    y += f32_cum_mul_scalar(&x0[N_vector], &x1[N_vector], N_remain);

    return y;
}
#endif

inline static 
float f32_cum_mul_auto(const float* x0, const float* x1, const int N) {
    #if defined(_DSP_AVX2)
    return f32_cum_mul_avx2(x0, x1, N);
    #elif defined(_DSP_SSSE3)
    return f32_cum_mul_ssse3(x0, x1, N);
    #else
    return f32_cum_mul_scalar(x0, x1, N);
    #endif
}