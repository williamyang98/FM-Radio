#pragma once
#include <cmath>
#include <assert.h>
#include <complex>

// NOTE: Assumes arrays are aligned
// Multiple a reference signal element-wise with a vector of complex floats
// The reference signal is given as:
// - Vector of floats representing time
// - Float representing multiple of time vector
// - Float representing offset of time vector

inline static 
void apply_harmonic_pll_scalar(
    const float* dt, const std::complex<float>* x, std::complex<float>* y, const int N,
    const float harmonic, const float offset) 
{
    for (int i = 0; i < N; i++) {
        const float dt_0 = dt[i];
        const float dt_c = dt_0*harmonic + offset;
        const auto pll = std::complex<float>(std::cos(dt_c), std::sin(dt_c));
        y[i] = x[i] * pll;
    }
}

// TODO: Modify code to support ARM platforms like Raspberry PI using NEON
#include <immintrin.h>
#include "data_packing.h"
#include "c32_mul.h"

#if !defined(_MSC_VER)
#pragma message("Compiling Harmonic PLL with external Intel SVML library")

#if defined(_DSP_SSSE3)
#define SSE_MATHFUN_WITH_CODE
#include "sse_mathfun.h"
#define _mm_cos_ps(x) cos_ps(x)
#endif

#if defined(_DSP_AVX2)
#include "avx_mathfun.h"
#define _mm256_cos_ps(x) cos256_ps(x)
#endif
#endif

#if defined(_DSP_SSSE3)
inline static
void apply_harmonic_pll_ssse3(
    const float* dt, const std::complex<float>* x, std::complex<float>* y, const int N,
    const float harmonic, const float offset) 
{
    // 128bits = 16bytes = 2*8bytes
    constexpr int K = 2;
    const int M = N/K;

    // Generate constants
    const __m128 harmonic_vec = _mm_set1_ps(harmonic);
    cpx128_t offset_vec;
    constexpr float PI = 3.14159265358979323846f;
    for (int i = 0; i < K; i++) {
        offset_vec.c32[i] = { 0.0f + offset, -PI/2.0f + offset };
    }

    // Constants for generating our dt vector
    // [3 2 1 0] -> [3 3 2 2]
    const uint8_t PERMUTE_UPPER = 0b11111010;
    // [3 2 1 0] -> [1 1 0 0]
    const uint8_t PERMUTE_LOWER = 0b01010000;

    for (int i = 0; i < M; i++) {
        // [c0 c1]
        __m128 b0 = _mm_load_ps(reinterpret_cast<const float*>(&x[i*K]));

        // [t0 t0 t1 t1]
        cpx128_t dt_vec; 
        dt_vec.f32[0] = dt[i*K+0];
        dt_vec.f32[1] = dt[i*K+0];
        dt_vec.f32[2] = dt[i*K+1];
        dt_vec.f32[3] = dt[i*K+1];

        // Apply harmonic
        dt_vec.ps = _mm_mul_ps(dt_vec.ps, harmonic_vec);
        // Apply offset and phase split for cos+jsin
        dt_vec.ps = _mm_add_ps(dt_vec.ps, offset_vec.ps);

        __m128 b1 = _mm_cos_ps(dt_vec.ps);
        __m128 res = c32_mul_ssse3(b0, b1);
        _mm_store_ps(reinterpret_cast<float*>(&y[i*K]), res);
    }

    const int N_vector = M*K;
    const int N_remain = N-N_vector;
    apply_harmonic_pll_scalar(
        &dt[N_vector], &x[N_vector], &y[N_vector], N_remain,
        harmonic, offset);
}
#endif

#if defined(_DSP_AVX2)
inline static 
void apply_harmonic_pll_avx2(
    const float* dt, const std::complex<float>* x, std::complex<float>* y, const int N,
    const float harmonic, const float offset) 
{
    // 256bits = 32bytes = 4*8bytes
    constexpr int K = 4;
    const int M = N/K;

    // Generate constants
    const __m256 harmonic_vec = _mm256_set1_ps(harmonic);
    cpx256_t offset_vec;
    constexpr float PI = 3.14159265358979323846f;
    for (int i = 0; i < K; i++) {
        offset_vec.c32[i] = { 0.0f + offset, -PI/2.0f + offset };
    }

    // Constants for generating our dt vector
    // [3 2 1 0] -> [3 3 2 2]
    const uint8_t PERMUTE_UPPER = 0b11111010;
    // [3 2 1 0] -> [1 1 0 0]
    const uint8_t PERMUTE_LOWER = 0b01010000;

    for (int i = 0; i < M; i++) {
        // [c0 c1 c2 c3]
        __m256 b0 = _mm256_load_ps(reinterpret_cast<const float*>(&x[i*K]));

        // [t0 t1 t2 t3]
        __m128 dt_sub_vec = _mm_load_ps(&dt[i*K]);
        // [t0 t0 t1 t1 t2 t2 t3 t3]
        cpx256_t dt_vec; 
        dt_vec.m128[0] = _mm_permute_ps(dt_sub_vec, PERMUTE_LOWER);
        dt_vec.m128[1] = _mm_permute_ps(dt_sub_vec, PERMUTE_UPPER);

        // Apply harmonic
        dt_vec.ps = _mm256_mul_ps(dt_vec.ps, harmonic_vec);
        // Apply offset and phase split for cos+jsin
        dt_vec.ps = _mm256_add_ps(dt_vec.ps, offset_vec.ps);

        __m256 b1 = _mm256_cos_ps(dt_vec.ps);
        __m256 res = c32_mul_avx2(b0, b1);
        _mm256_store_ps(reinterpret_cast<float*>(&y[i*K]), res);
    }

    const int N_vector = M*K;
    const int N_remain = N-N_vector;
    apply_harmonic_pll_scalar(
        &dt[N_vector], &x[N_vector], &y[N_vector], N_remain,
        harmonic, offset);
}
#endif

inline static
void apply_harmonic_pll_auto(
    const float* dt, const std::complex<float>* x, std::complex<float>* y, const int N,
    const float harmonic, const float offset) 
{
    #if defined(_DSP_AVX2)
    apply_harmonic_pll_avx2(dt, x, y, N, harmonic, offset);
    #elif defined(_DSP_SSSE3)
    apply_harmonic_pll_ssse3(dt, x, y, N, harmonic, offset);
    #else
    apply_harmonic_pll_scalar(dt, x, y, N, harmonic, offset);
    #endif
}