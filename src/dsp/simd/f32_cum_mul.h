#pragma once
#include <assert.h>

// NOTE: Assumes arrays are aligned
// Multiply and accumulate vector of floats with another vector of floats



#pragma once

#include "./detect_architecture.h"
#include "./dsp_config.h"
#include "./scalar/f32_cum_mul.h"
#if defined(__ARCH_X86__)
#include "./x86/f32_cum_mul.h"
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