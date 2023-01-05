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
#include "f32_cum_sum.h"

#if defined(_DSP_SSSE3)
static inline
float f32_cum_mul_ssse3(const float* x0, const float* x1, const int N)
{
    float y = 0;

    // 128bits = 16bytes = 4*4bytes
    const int K = 4;
    const int M = N/K;

    cpx128_t v_sum;
    v_sum.ps = _mm_set1_ps(0.0f);

    for (int i = 0; i < M; i++) {
        __m128 a0 = _mm_load_ps(&x0[i*K]);
        __m128 a1 = _mm_load_ps(&x1[i*K]);

        // multiply accumulate
        #if !defined(_DSP_FMA)
        v_sum.ps = _mm_add_ps(_mm_mul_ps(a0, a1), v_sum.ps);
        #else
        v_sum.ps = _mm_fmadd_ps(a0, a1, v_sum.ps);
        #endif
    }

    y += f32_cum_sum_ssse3(v_sum);

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

    cpx256_t v_sum;
    v_sum.ps = _mm256_set1_ps(0.0f);

    for (int i = 0; i < M; i++) {
        __m256 a0 = _mm256_load_ps(&x0[i*K]);
        __m256 a1 = _mm256_load_ps(&x1[i*K]);

        // multiply accumulate
        #if !defined(_DSP_FMA)
        v_sum.ps = _mm256_add_ps(_mm256_mul_ps(a0, a1), v_sum.ps);
        #else
        v_sum.ps = _mm256_fmadd_ps(a0, a1, v_sum.ps);
        #endif
    }

    y += f32_cum_sum_avx2(v_sum);

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