#define _USE_MATH_DEFINES
#include <cmath>

#include "hilbert_fir_filter.h"

#define _min(A,B) (A > B) ? B : A
#define _max(A,B) (A > B) ? A : B

Hilbert_FIR_Filter::Hilbert_FIR_Filter(const int _K)
: K(_K), 
  xn(_K), b(_K)
{
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
            // NOTE: time reverse coefficients here
            const float _n = -(float)n;  
            const float _a = std::sin(_n*PI/2.0f);
            b[i] = 2.0f/(_n*PI) * _a * _a;
        }
    }

    for (int i = 0; i < K; i++) {
        xn[i] = 0;
    }
}

void Hilbert_FIR_Filter::process(const float* x, std::complex<float>* y, const int N)
{
    const int M = (K-1)/2;

    // continue from previous block
    const int M0 = _min(K-1, N);
    for (int i = 0; i < M0; i++) {
        push_value(x[i]);
        const float imag = apply_filter(xn.data());
        y[i] = { xn[M], imag };
    }

    // inplace math
    for (int i = M0, j = 0; i < N; i++, j++) {
        const float imag = apply_filter(&x[j]);
        y[i] = { x[j+M], imag };
    }

    // push end of buffer
    const int M1 = _max(N-K, M0);
    push_values(&x[M1], N-M1);
}

void Hilbert_FIR_Filter::process(const float* x, float* y, const int N)
{
    const int M = (K-1)/2;

    // continue from previous block
    const int M0 = _min(K-1, N);
    for (int i = 0; i < M0; i++) {
        push_value(x[i]);
        y[i] = apply_filter(xn.data());
    }

    // inplace math
    for (int i = M0, j = 0; i < N; i++, j++) {
        y[i] = apply_filter(&x[j]);
    }

    // push end of buffer
    const int M1 = _max(N-K, M0);
    push_values(&x[M1], N-M1);
}