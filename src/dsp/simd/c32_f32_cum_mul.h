#pragma once

#include <complex>
#include "./detect_architecture.h"
#include "./dsp_config.h"
#include "./scalar/c32_f32_cum_mul.h"
#if defined(__ARCH_X86__)
#include "./x86/c32_f32_cum_mul.h"
#endif

inline static 
std::complex<float> c32_f32_cum_mul_auto(const std::complex<float>* x0, const float* x1, const int N) {
    #if defined(_DSP_AVX2)
    return c32_f32_cum_mul_avx2(x0, x1, N);
    #elif defined(_DSP_SSSE3)
    return c32_f32_cum_mul_ssse3(x0, x1, N);
    #else
    return c32_f32_cum_mul_scalar(x0, x1, N);
    #endif
}