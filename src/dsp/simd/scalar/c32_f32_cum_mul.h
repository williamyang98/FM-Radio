#pragma once
#include <complex>

static inline
std::complex<float> c32_f32_cum_mul_scalar(const std::complex<float>* x0, const float* x1, const int N) {
    auto y = std::complex<float>(0,0);
    for (int i = 0; i < N; i++) {
        y += x0[i] * x1[i];
    }
    return y;
}