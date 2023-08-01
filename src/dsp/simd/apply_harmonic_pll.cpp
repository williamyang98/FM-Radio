#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <assert.h>
#include <stdalign.h>

#include "./apply_harmonic_pll.h"
#include "./detect_architecture.h"
#include "./simd_flags.h"

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

#if defined(__ARCH_X86__)

#if defined(_MSC_VER) && !defined(__clang__)
#define IS_INTEL_SVML_PRESENT 1
#endif

#if defined(__SSE3__)
#include <immintrin.h>
#include "./x86/c32_mul.h"

#if !IS_INTEL_SVML_PRESENT
#define USE_SSE_AUTO
#define SSE_MATHFUN_WITH_CODE
#include "./x86/sse_mathfun.h"
#define _mm_cos_ps(x) cos_ps(x)
#endif

void apply_harmonic_pll_sse3(
    const float* dt, const std::complex<float>* x, std::complex<float>* y, const int N,
    const float harmonic, const float offset) 
{
    // 128bits = 16bytes = 2*8bytes
    constexpr int K = 2;
    const int M = N/K;
    const int N_vector = M*K;
    const int N_remain = N-N_vector;

    // Generate constants
    const __m128 harmonic_vec = _mm_set1_ps(harmonic);
    alignas(16) std::complex<float> _offset_vec[K];
    for (int i = 0; i < K; i++) {
        _offset_vec[i] = { 0.0f + offset, -float(M_PI)/2.0f + offset };
    }
    __m128 offset_vec = _mm_load_ps(reinterpret_cast<const float *>(_offset_vec));

    // Constants for generating our dt vector
    // [3 2 1 0] -> [3 3 2 2]
    const uint8_t PERMUTE_UPPER = 0b11111010;
    // [3 2 1 0] -> [1 1 0 0]
    const uint8_t PERMUTE_LOWER = 0b01010000;

    for (int i = 0; i < N_vector; i+=K) {
        // [c0 c1]
        __m128 b0 = _mm_loadu_ps(reinterpret_cast<const float*>(&x[i]));

        // [t0 t0 t1 t1]
        __m128 dt_vec = _mm_set_ps(
            dt[i+1], dt[i+1],
            dt[i+0], dt[i+0]
        );

        // Apply harmonic
        dt_vec = _mm_mul_ps(dt_vec, harmonic_vec);
        // Apply offset and phase split for cos+jsin
        dt_vec = _mm_add_ps(dt_vec, offset_vec);

        __m128 b1 = _mm_cos_ps(dt_vec);
        __m128 res = c32_mul_sse3(b0, b1);
        _mm_storeu_ps(reinterpret_cast<float*>(&y[i]), res);
    }

    apply_harmonic_pll_scalar(
        &dt[N_vector], &x[N_vector], &y[N_vector], N_remain,
        harmonic, offset);
}
#endif

#if defined(__AVX__)
#include <immintrin.h>
#include "./x86/c32_mul.h"

#if !IS_INTEL_SVML_PRESENT
#include "./x86/avx_mathfun.h"
#define _mm256_cos_ps(x) cos256_ps(x)
#endif

void apply_harmonic_pll_avx(
    const float* dt, const std::complex<float>* x, std::complex<float>* y, const int N,
    const float harmonic, const float offset) 
{
    // 256bits = 32bytes = 4*8bytes
    constexpr int K = 4;
    const int M = N/K;
    const int N_vector = M*K;
    const int N_remain = N-N_vector;

    // Generate constants
    const __m256 harmonic_vec = _mm256_set1_ps(harmonic);
    alignas(32) std::complex<float> _offset_vec[K];
    for (int i = 0; i < K; i++) {
        _offset_vec[i] = { 0.0f + offset, -float(M_PI)/2.0f + offset };
    }
    __m256 offset_vec = _mm256_load_ps(reinterpret_cast<const float*>(_offset_vec));

    // Constants for generating our dt vector
    // [3 2 1 0] -> [3 3 2 2]
    const uint8_t PERMUTE_UPPER = 0b11111010;
    // [3 2 1 0] -> [1 1 0 0]
    const uint8_t PERMUTE_LOWER = 0b01010000;

    for (int i = 0; i < N_vector; i+=K) {
        // [c0 c1 c2 c3]
        __m256 b0 = _mm256_loadu_ps(reinterpret_cast<const float*>(&x[i]));

        // [t0 t1 t2 t3]
        __m128 dt_sub_vec = _mm_loadu_ps(&dt[i]);
        // [t0 t0 t1 t1 t2 t2 t3 t3]
        __m256 dt_vec = _mm256_set_m128(
            _mm_permute_ps(dt_sub_vec, PERMUTE_UPPER),
            _mm_permute_ps(dt_sub_vec, PERMUTE_LOWER)
        );

        // Apply harmonic
        dt_vec = _mm256_mul_ps(dt_vec, harmonic_vec);
        // Apply offset and phase split for cos+jsin
        dt_vec = _mm256_add_ps(dt_vec, offset_vec);

        __m256 b1 = _mm256_cos_ps(dt_vec);
        __m256 res = c32_mul_avx(b0, b1);
        _mm256_storeu_ps(reinterpret_cast<float*>(&y[i]), res);
    }

    apply_harmonic_pll_scalar(
        &dt[N_vector], &x[N_vector], &y[N_vector], N_remain,
        harmonic, offset);
}
#endif

#endif

void apply_harmonic_pll_auto(
    const float* dt, const std::complex<float>* x, std::complex<float>* y, const int N,
    const float harmonic, const float offset) 
{
    #if defined(__ARCH_X86__)
        #if defined(__AVX__)
            apply_harmonic_pll_avx(dt, x, y, N, harmonic, offset);
        #elif defined(__SSE3__)
            apply_harmonic_pll_sse3(dt, x, y, N, harmonic, offset);
        #else
            apply_harmonic_pll_scalar(dt, x, y, N, harmonic, offset);
        #endif
    #else
        apply_harmonic_pll_scalar(dt, x, y, N, harmonic, offset);
    #endif
}
