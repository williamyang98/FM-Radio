#pragma once

#include <complex>

// Implement the Hilbert transform using an FIR filter
class Hilbert_FIR_Filter 
{
public:
    float* xn;
    float* b;
    const int K;
    int Ki;
public:
    Hilbert_FIR_Filter(const int _K);
    ~Hilbert_FIR_Filter();
    // Generate the quadrature component
    // Real component is delayed by floor(K/2) samples to match delay of hilbert filter
    // Imag component contains the Hilbert transform delayed by floor(K/2) samples
    void process(const float* x, std::complex<float>* y, const int N);
    // Generate the shifted real component
    void process(const float* x, float* y, const int N);
};
