#pragma once

#include "utility/aligned_vector.h"
#include <complex>

// Implement the Hilbert transform using an FIR filter
class Hilbert_FIR_Filter 
{
public:
    AlignedVector<float> xn;
    AlignedVector<float> b;
    const int K;
public:
    Hilbert_FIR_Filter(const int _K);
    // Generate the quadrature component
    // Real component is delayed by floor(K/2) samples to match delay of hilbert filter
    // Imag component contains the Hilbert transform delayed by floor(K/2) samples
    void process(const float* x, std::complex<float>* y, const int N);
    // Generate the shifted real component
    void process(const float* x, float* y, const int N);
private:
    void push_value(float x) {
        for (int i = 0; i < (K-1); i++) {
            xn[i] = xn[i+1];
        }
        xn[K-1] = x;
    }

    void push_values(const float* x, const int N) {
        const int M = K-N;
        for (int i = 0; i < M; i++) {
            xn[i] = xn[i+N];
        }
        for (int i = M, j = 0; i < K; i++, j++) {
            xn[i] = x[j];
        }
    }

    float apply_filter(const float* x) {
        float y = 0; 
        for (int i = 0; i < K; i++) {
            y += x[i] * b[i];
        }
        return y;
    }
};
