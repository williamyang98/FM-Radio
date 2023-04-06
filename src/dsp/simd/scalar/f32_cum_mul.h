#pragma once

static inline
float f32_cum_mul_scalar(const float* x0, const float* x1, const int N) {
    float y = 0;
    for (int i = 0; i < N; i++) {
        y += x0[i] * x1[i];
    }
    return y;
}