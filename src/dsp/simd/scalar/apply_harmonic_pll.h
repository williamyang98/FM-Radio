#pragma once
#include <cmath>
#include <assert.h>
#include <complex>

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