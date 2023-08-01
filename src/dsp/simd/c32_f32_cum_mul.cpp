#include "./c32_f32_cum_mul.h"
#include "./detect_architecture.h"
#include "./simd_flags.h"

std::complex<float> c32_f32_cum_mul_scalar(const std::complex<float>* x0, const float* x1, const int N) {
    auto y = std::complex<float>(0,0);
    for (int i = 0; i < N; i++) {
        y += x0[i] * x1[i];
    }
    return y;
}

#if defined(__ARCH_X86__)

#if defined(__SSE2__)
#include <immintrin.h>
#include "./x86/c32_mul.h"
#include "./x86/c32_cum_sum.h"

std::complex<float> c32_f32_cum_mul_sse2(const std::complex<float>* x0, const float* x1, const int N)
{
    // 128bits = 16bytes = 4*4bytes
    constexpr int K = 4;
    const int M = N/K;
    const int N_vector = M*K;
    const int N_remain = N-N_vector;

    // [3 2 1 0] -> [3 3 2 2]
    const uint8_t PERMUTE_UPPER = 0b11111010;
    // [3 2 1 0] -> [1 1 0 0]
    const uint8_t PERMUTE_LOWER = 0b01010000;

    __m128 v_sum = _mm_set1_ps(0.0f);

    for (int i = 0; i < N_vector; i+=K) {
        // [c0 c1]
        __m128 a0 = _mm_loadu_ps(reinterpret_cast<const float*>(&x0[i]));
        // [c2 c3]
        __m128 a1 = _mm_loadu_ps(reinterpret_cast<const float*>(&x0[i + K/2]));

        // [a0 a1 a2 a3]
        __m128 b0 = _mm_loadu_ps(&x1[i]);
        // [a2 a2 a3 a3]
        __m128 b1 = _mm_shuffle_ps(b0, b0, PERMUTE_UPPER);
        // [a0 a0 a1 a1]
        b0 = _mm_shuffle_ps(b0, b0, PERMUTE_LOWER);

        // multiply accumulate
        #if !defined(__FMA__)
        v_sum = _mm_add_ps(_mm_mul_ps(a0, b0), v_sum);
        v_sum = _mm_add_ps(_mm_mul_ps(a1, b1), v_sum);
        #else
        v_sum = _mm_fmadd_ps(a0, b0, v_sum);
        v_sum = _mm_fmadd_ps(a1, b1, v_sum);
        #endif
    }

    auto y = std::complex<float>(0,0);
    y += c32_cum_sum_sse2(v_sum);
    y += c32_f32_cum_mul_scalar(&x0[N_vector], &x1[N_vector], N_remain);
    return y;
}
#endif

#if defined(__AVX__)
#include <immintrin.h>
#include "./x86/c32_mul.h"
#include "./x86/c32_cum_sum.h"

std::complex<float> c32_f32_cum_mul_avx(const std::complex<float>* x0, const float* x1, const int N)
{
    // 256bits = 32bytes = 4*8bytes
    constexpr int K = 4;
    const int M = N/K;
    const int N_vector = M*K;
    const int N_remain = N-N_vector;

    // [3 2 1 0] -> [3 3 2 2]
    const uint8_t PERMUTE_UPPER = 0b11111010;
    // [3 2 1 0] -> [1 1 0 0]
    const uint8_t PERMUTE_LOWER = 0b01010000;

    __m256 v_sum = _mm256_set1_ps(0.0f);

    for (int i = 0; i < N_vector; i+=K) {
        // [c0 c1 c2 c3]
        __m256 a0 = _mm256_loadu_ps(reinterpret_cast<const float*>(&x0[i]));

        // [a0 a1 a2 a3]
        __m128 b0 = _mm_loadu_ps(&x1[i]);
        // [a2 a2 a3 a3]
        __m128 b1 = _mm_permute_ps(b0, PERMUTE_LOWER);
        // [a0 a0 a1 a1]
        b0 = _mm_permute_ps(b0, PERMUTE_UPPER);

        // [a0 a0 a1 a1 a2 a2 a3 a3]
        __m256 a1 = _mm256_set_m128(b0, b1);

        // multiply accumulate
        #if !defined(__FMA__)
        v_sum = _mm256_add_ps(_mm256_mul_ps(a0, a1), v_sum);
        #else
        v_sum = _mm256_fmadd_ps(a0, a1, v_sum);
        #endif
    }

    auto y = std::complex<float>(0,0);
    y += c32_cum_sum_avx(v_sum);
    y += c32_f32_cum_mul_scalar(&x0[N_vector], &x1[N_vector], N_remain);
    return y;
}
#endif

#endif

std::complex<float> c32_f32_cum_mul_auto(const std::complex<float>* x0, const float* x1, const int N) {
    #if defined(__ARCH_X86__)
        #if defined(__AVX__)
            return c32_f32_cum_mul_avx(x0, x1, N);
        #elif defined(__SSE2__)
            return c32_f32_cum_mul_sse2(x0, x1, N);
        #else
            return c32_f32_cum_mul_scalar(x0, x1, N);
        #endif
    #else
        return c32_f32_cum_mul_scalar(x0, x1, N);
    #endif
}
