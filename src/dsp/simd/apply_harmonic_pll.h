#pragma once

#include <complex>
void apply_harmonic_pll_auto(
    const float* dt, const std::complex<float>* x, std::complex<float>* y, const int N,
    const float harmonic, const float offset);