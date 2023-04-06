#pragma once

#include <complex>
#include "./detect_architecture.h"
#include "./dsp_config.h"
#include "./scalar/apply_harmonic_pll.h"
#if defined(__ARCH_X86__)
#include "./x86/apply_harmonic_pll.h"
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