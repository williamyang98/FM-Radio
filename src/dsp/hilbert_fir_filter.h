#pragma once

#include "utility/aligned_vector.h"
#include <complex>
#include <cmath>

#define _min(A,B) (A > B) ? B : A
#define _max(A,B) (A > B) ? A : B

// Implement the Hilbert transform using an FIR filter
template <typename T>
class Hilbert_FIR_Filter 
{
public:
    const int K;
    AlignedVector<float> b;
    AlignedVector<T> xn;
    AlignedVector<T> tmp;
public:
    Hilbert_FIR_Filter(const int _K)
    : K(_K), 
      b(_K), xn(_K), tmp(_K)
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
    // Generate the quadrature component
    // Real component is delayed by floor(K/2) samples to match delay of hilbert filter
    // Imag component contains the Hilbert transform delayed by floor(K/2) samples
    void process(const T* x, std::complex<T>* y, const int N) {
        const int M = (K-1)/2;
        const int M0 = _min(K-1, N);
        const int M1 = _max(N-K, M0);

        // save the tail which is overridden
        for (int i = M1, j = 0; i < N; i++, j++) {
            tmp[j] = x[i];
        }

        // inplace math
        // NOTE: We do this first and reverse so we don't destroy the first K samples
        //       This can occur if x and y are the same buffer
        for (int i = M0, j = 0; i < N; i++, j++) {
            const T imag = apply_filter(&x[j]);
            y[i] = { x[j+M], imag };
        }

        // continue from previous block
        for (int i = 0; i < M0; i++) {
            push_value(x[i]);
            const T imag = apply_filter(xn.data());
            y[i] = { xn[M], imag };
        }

        // push end of buffer
        push_values(tmp.data(), N-M1);
    }
    // Generate the shifted real component
    void process(const T* x, T* y, const int N) {
        const int M0 = _min(K-1, N);
        const int M1 = _max(N-K, M0);

        // save the tail which is overridden
        for (int i = M1, j = 0; i < N; i++, j++) {
            tmp[j] = x[i];
        }

        // inplace math
        // NOTE: We do this first and reverse so we don't destroy the first K samples
        //       This can occur if x and y are the same buffer
        for (int i = (N-1), j = (N-1)-(K-1); i >= M0; i--, j--) {
            y[i] = apply_filter(&x[j]);
        }

        // continue from previous block
        for (int i = 0; i < M0; i++) {
            push_value(x[i]);
            y[i] = apply_filter(xn.data());
        }

        // push end of buffer
        push_values(tmp.data(), N-M1);
    }
private:
    void push_value(T x) {
        for (int i = 0; i < (K-1); i++) {
            xn[i] = xn[i+1];
        }
        xn[K-1] = x;
    }

    void push_values(const T* x, const int N) {
        const int M = K-N;
        for (int i = 0; i < M; i++) {
            xn[i] = xn[i+N];
        }
        for (int i = M, j = 0; i < K; i++, j++) {
            xn[i] = x[j];
        }
    }

    T apply_filter(const T* x) {
        T y;
        y = 0; 
        for (int i = 0; i < K; i++) {
            y += x[i] * b[i];
        }
        return y;
    }
};

#undef _min
#undef _max

#include "simd/f32_cum_mul.h"
float Hilbert_FIR_Filter<float>::apply_filter(const float* x) {
    return f32_cum_mul_auto(x, b.data(), K);
}