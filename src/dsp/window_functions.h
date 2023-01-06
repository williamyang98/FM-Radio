#pragma once
#include <cmath>

// Window functions
// For a buffer with N samples
// -N/2 <= i <= +N/2
// M = N-1
// x = 2*pi*i/M

static inline 
float window_hamming(const float x) {
    return 0.53836f - 0.46164f*std::cos(x);
}

static inline
float window_hann(const float x) {
    const float a = std::sin(x/2.0f);
    return a*a;
}

static inline
float window_blackman(const float x) {
    const float a0 = 0.42659f;
    const float a1 = 0.49656f;
    const float a2 = 0.076849f;
    return a0 - a1*std::cos(x) + a2*std::cos(2.0f*x);
}

static inline
float window_blackman_harris(const float x) {
    const float a0 = 0.35875f;
    const float a1 = 0.48829f;
    const float a2 = 0.14128f;
    const float a3 = 0.01168f;
    return a0 - a1*std::cos(x) + a2*std::cos(2.0f*x) - a3*std::cos(3.0f*x);
}