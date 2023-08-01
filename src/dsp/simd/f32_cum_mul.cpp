#include "./f32_cum_mul.h"
#include "./detect_architecture.h"
#include "./simd_flags.h"

float f32_cum_mul_scalar(const float* x0, const float* x1, const int N) {
    float y = 0;
    for (int i = 0; i < N; i++) {
        y += x0[i] * x1[i];
    }
    return y;
}

#if defined(__ARCH_X86__)

#if defined(__SSE2__)
#include <immintrin.h>
#include "./x86/f32_cum_sum.h"

float f32_cum_mul_sse2(const float* x0, const float* x1, const int N)
{
    // 128bits = 16bytes = 4*4bytes
    const int K = 4;
    const int M = N/K;
    const int N_vector = M*K;
    const int N_remain = N-N_vector;

    __m128 v_sum = _mm_set1_ps(0.0f);

    for (int i = 0; i < N_vector; i+=K) {
        __m128 a0 = _mm_loadu_ps(&x0[i]);
        __m128 a1 = _mm_loadu_ps(&x1[i]);

        // multiply accumulate
        #if !defined(__FMA__)
        v_sum = _mm_add_ps(_mm_mul_ps(a0, a1), v_sum);
        #else
        v_sum = _mm_fmadd_ps(a0, a1, v_sum);
        #endif
    }

    float y = 0;
    y += f32_cum_sum_sse2(v_sum);
    y += f32_cum_mul_scalar(&x0[N_vector], &x1[N_vector], N_remain);
    return y;
}
#endif

#if defined(__AVX__)
#include <immintrin.h>
#include "./x86/f32_cum_sum.h"

float f32_cum_mul_avx(const float* x0, const float* x1, const int N)
{
    // 256bits = 32bytes = 8*4bytes
    const int K = 8;
    const int M = N/K;
    const int N_vector = M*K;
    const int N_remain = N-N_vector;

    __m256 v_sum = _mm256_set1_ps(0.0f);

    for (int i = 0; i < N_vector; i+=K) {
        __m256 a0 = _mm256_loadu_ps(&x0[i]);
        __m256 a1 = _mm256_loadu_ps(&x1[i]);

        // multiply accumulate
        #if !defined(__FMA__)
        v_sum = _mm256_add_ps(_mm256_mul_ps(a0, a1), v_sum);
        #else
        v_sum = _mm256_fmadd_ps(a0, a1, v_sum);
        #endif
    }

    float y = 0;
    y += f32_cum_sum_avx(v_sum);
    y += f32_cum_mul_scalar(&x0[N_vector], &x1[N_vector], N_remain);
    return y;
}
#endif

#endif

float f32_cum_mul_auto(const float* x0, const float* x1, const int N) {
    #if defined(__ARCH_X86__)
        #if defined(__AVX__)
            return f32_cum_mul_avx(x0, x1, N);
        #elif defined(__SSE2__)
            return f32_cum_mul_sse2(x0, x1, N);
        #else
            return f32_cum_mul_scalar(x0, x1, N);
        #endif
    #else
        return f32_cum_mul_scalar(x0, x1, N);
    #endif
}
