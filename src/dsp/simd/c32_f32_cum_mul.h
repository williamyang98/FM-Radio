#pragma once
#include <assert.h>
#include <complex>

// NOTE: Assumes arrays are aligned
// Multiply and accumulate vector of complex floats with vector of floats

static inline
std::complex<float> c32_f32_cum_mul_scalar(const std::complex<float>* x0, const float* x1, const int N) {
    auto y = std::complex<float>(0,0);
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
std::complex<float> c32_f32_cum_mul_ssse3(const std::complex<float>* x0, const float* x1, const int N)
{
    auto y = std::complex<float>(0,0);

    // 128bits = 16bytes = 4*4bytes
    constexpr int K = 4;
    const int M = N/K;

    // [3 2 1 0] -> [3 3 2 2]
    const uint8_t PERMUTE_UPPER = 0b11111010;
    // [3 2 1 0] -> [1 1 0 0]
    const uint8_t PERMUTE_LOWER = 0b01010000;

    for (int i = 0; i < M; i++) {
        // [c0 c1]
        __m128 a0 = _mm_load_ps(reinterpret_cast<const float*>(&x0[i*K]));
        // [c2 c3]
        __m128 a1 = _mm_load_ps(reinterpret_cast<const float*>(&x0[i*K + K/2]));

        // [a0 a1 a2 a3]
        __m128 b0 = _mm_load_ps(&x1[i*K]);
        // [a2 a2 a3 a3]
        __m128 b1 = _mm_shuffle_ps(b0, b0, PERMUTE_UPPER);
        // [a0 a0 a1 a1]
        b0 = _mm_shuffle_ps(b0, b0, PERMUTE_LOWER);

        // multiply vector
        cpx128_t c0, c1;
        c0.ps = _mm_mul_ps(a0, b0);
        c1.ps = _mm_mul_ps(a1, b1);

        // Perform vectorised cumulative sum
        // Shift half of vector and add. Repeat until we get the final sum
        cpx128_t d0, d1;
        // [c0+c2 c1+c3]
        d0.ps = _mm_add_ps(c0.ps, c1.ps);
        // [c1+c3 0]
        d1.i = _mm_srli_si128(d0.i, 1*sizeof(std::complex<float>));
        // [c0+c1+c2+c3]
        d0.ps = _mm_add_ps(d0.ps, d1.ps);

        y += d0.c32[0];
    }

    const int N_vector = M*K;
    const int N_remain = N-N_vector;
    y += c32_f32_cum_mul_scalar(&x0[N_vector], &x1[N_vector], N_remain);

    return y;
}
#endif

#if defined(_DSP_AVX2)
static inline
std::complex<float> c32_f32_cum_mul_avx2(const std::complex<float>* x0, const float* x1, const int N)
{
    auto y = std::complex<float>(0,0);

    // 256bits = 32bytes = 4*8bytes
    constexpr int K = 4;
    const int M = N/K;

    // [3 2 1 0] -> [3 3 2 2]
    const uint8_t PERMUTE_UPPER = 0b11111010;
    // [3 2 1 0] -> [1 1 0 0]
    const uint8_t PERMUTE_LOWER = 0b01010000;

    for (int i = 0; i < M; i++) {
        // [c0 c1 c2 c3]
        __m256 a0 = _mm256_load_ps(reinterpret_cast<const float*>(&x0[i*K]));

        // [a0 a1 a2 a3]
        __m128 b0 = _mm_load_ps(&x1[i*K]);
        // [a2 a2 a3 a3]
        __m128 b1 = _mm_permute_ps(b0, PERMUTE_LOWER);
        // [a0 a0 a1 a1]
        b0 = _mm_permute_ps(b0, PERMUTE_UPPER);

        // [a0 a0 a1 a1 a2 a2 a3 a3]
        __m256 a1 = _mm256_set_m128(b0, b1);

        // multiply vector
        cpx256_t c0;
        c0.ps = _mm256_mul_ps(a0, a1);

        // Perform vectorised cumulative sum
        // Shift half of vector and add. Repeat until we get the final sum
        cpx128_t d0, d1;
        // [c0+c2 c1+c3]
        d0.ps = _mm_add_ps(c0.m128[0], c0.m128[1]);
        // [c1+c3 0]
        d1.i = _mm_srli_si128(d0.i, 1*sizeof(std::complex<float>));
        // [c0+c1+c2+c3 0]
        d0.ps = _mm_add_ps(d0.ps, d1.ps);

        y += d0.c32[0];
    }

    const int N_vector = M*K;
    const int N_remain = N-N_vector;
    y += c32_f32_cum_mul_scalar(&x0[N_vector], &x1[N_vector], N_remain);

    return y;
}
#endif

inline static 
std::complex<float> c32_f32_cum_mul_auto(const std::complex<float>* x0, const float* x1, const int N) {
    #if defined(_DSP_AVX2)
    // return c32_f32_cum_mul_avx2(x0, x1, N);
    // FIXME: AVX2 performs worse than SSE3 for unknown reasons
    return c32_f32_cum_mul_ssse3(x0, x1, N);
    #elif defined(_DSP_SSSE3)
    return c32_f32_cum_mul_ssse3(x0, x1, N);
    #else
    return c32_f32_cum_mul_scalar(x0, x1, N);
    #endif
}
