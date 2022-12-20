#define _USE_MATH_DEFINES
#include <cmath>

#include "hilbert_fir_filter.h"

Hilbert_FIR_Filter::Hilbert_FIR_Filter(const int _K)
: K(_K)
{
    xn = new float[K]{0};
    b = new float[K];
    Ki = 0;

    // Generate the Hilbert filter
    // Given by the non-causal equation 
    // n != 0 : 2/(n*pi) * sin(n*pi/2)^2 
    // n  = 0 : 0
    const int M = (K-1)/2;
    for (int i = 0; i < K; i++) {
        const int n = i-M;
        if (n == 0) {
            b[i] = 0.0f;
        } else {
            constexpr float PI = (float)M_PI;
            const float _n = (float)n;
            const float _a = std::sin(_n*PI/2.0f);
            b[i] = 2.0f/(_n*PI) * _a * _a;
        }
    }
}

Hilbert_FIR_Filter::~Hilbert_FIR_Filter() {
    delete [] b;
    delete [] xn;
}

void Hilbert_FIR_Filter::process(const float* x, std::complex<float>* y, const int N)
{
    const int M = (K-1)/2;
    for (int i = 0; i < N; i++) {
        xn[Ki] = x[i];
        float imag = 0.0f;
        for (int j = 0; j < K; j++) {
            imag += xn[(K+Ki-j)%K] * b[j];
        }

        const int j = (K+Ki-M)%K;
        y[i] = { xn[j], imag };
        Ki = (Ki+1)%K;
    }
}

void Hilbert_FIR_Filter::process(const float* x, float* y, const int N)
{
    for (int i = 0; i < N; i++) {
        xn[Ki] = x[i];
        y[i] = 0.0f;
        for (int j = 0; j < K; j++) {
            y[i] += xn[(K+Ki-j)%K] * b[j];
        }
        Ki = (Ki+1)%K;
    }
}