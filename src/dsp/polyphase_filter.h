#pragma once
#include <vector>
#include "utility/aligned_allocator.hpp"

#define _min(A,B) (A > B) ? B : A
#define _max(A,B) (A > B) ? A : B

template <typename T>
class PolyphaseDownsampler
{
private:
    const int M;
    const int K;
    const int NN;
    std::vector<float, AlignedAllocator<float>> b;
    std::vector<T, AlignedAllocator<T>> xn;
public:
    float* get_b() { return b.data(); }
    int get_K() const { return NN; }
public:
    // b = FIR filter with M*K coefficients
    // M = downsampling factor and total phases 
    // K = total coefficients per phase
    PolyphaseDownsampler(const int _M, const int _K) 
    : M(_M), K(_K), NN(_M*_K),
      b(NN), xn(NN)
     {
        for (int i = 0; i < NN; i++) {
            b[i] = 0;
            xn[i] = 0;
        }
    }

    // E.g. M = 3, K = 2, M*K = 6
    // x8 x7 x6 x5 x4 x3 x2 x1 x0
    // b5 b4 b3 b2 b1 b0          => y0
    //    b5 b4 b3 b2 b1 b0       => ... (downsampling ignores these)
    //       b5 b4 b3 b2 b1 b0    => ...
    //          b5 b4 b3 b2 b1 b0 => y1
    // N = produce N output samples
    void process(const T* x, T* y, const int N) {
        const int M0 = _min(K-1, N);

        // NOTE: When downsampling we don't expect x and y to be the same buffer
        // continue from previous block
        for (int i = 0, j = 0; i < M0; i++, j+=M) {
            push_values(&x[j], M);
            y[i] = apply_filter(xn.data());
        }


        // inplace math
        for (int i = M0, j = 0; i < N; i++, j+=M) {
            y[i] = apply_filter(&x[j]);
        }

        // push end of buffer
        const int M1 = _max(N-K, M0);
        {
            const int i = M1;
            const int j = i*M;
            push_values(&x[j], (N-i)*M);
        }
    }

private:
    inline
    void push_values(const T* x, const int N) {
        const int M = NN-N;
        for (int i = 0; i < M; i++) {
            xn[i] = xn[i+N];
        }
        for (int i = M, j = 0; i < NN; i++, j++) {
            xn[i] = x[j];
        }
    }

    inline
    T apply_filter(const T* x) {
        T y;
        y = 0;
        for (int i = 0; i < NN; i++) {
            y += x[i] * b[i];
        }
        return y;
    }
};

template <typename T>
class PolyphaseUpsampler
{
private:
    const int L;
    const int K;
    const int NN;
    std::vector<float, AlignedAllocator<float>> b;
    std::vector<T, AlignedAllocator<T>> xn;
public:
    // b = FIR filter coefficients of length L*K
    // L = upsampling factor and total phases
    // K = number of coefficients per phase
    PolyphaseUpsampler(const float *_b, const int _L, const int _K) 
    : L(_L), K(_K), NN(_L*_K),
      b(NN), xn(K)
    {
        // TODO: Determine if we can use filter designer without repacking coefficients
        for (int phase = 0; phase < L; phase++) {
            const int phase_c = (L-1)-phase;
            for (int i = 0; i < K; i++) {
                const int j0 = phase_c*K + i;
                const int j1 = phase + i*L;
                // repack the coefficients so that they are contiguous
                b[j0] = _b[(NN-1)-j1] * (float)L;
            }
        }

        for (int i = 0; i < K; i++) {
            xn[i] = 0;
        }
    }

    // E.g. L = 4, K = 2, L*K = 8
    //  0  0  0 x0  0  0  0 x1  0  0  0 x2 
    // b7 b6 b5 b4 b3 b2 b1 b0          => y0
    //    b7 b6 b5 b4 b3 b2 b1 b0       => y1
    //       b7 b6 b5 b4 b3 b2 b1 b0    => y2
    //          b7 b6 b5 b4 b3 b2 b1 b0 => y3
    // N = process N input samples
    void process(const T* x, T* y, const int N) {
        const int M0 = _min(K-1, N);

        // NOTE: When upsampling we don't expect x and y to be the same buffer
        // continue from previous block
        for (int i = 0; i < M0; i++) {
            push_value(x[i]);
            for (int phase = 0; phase < L; phase++) {
                const int yi = i*L + phase;
                y[yi] = apply_filter(xn.data(), phase);
            }
        }

        // inplace math
        for (int i = M0, j = 0; i < N; i++, j++) {
            for (int phase = 0; phase < L; phase++) {
                const int yi = i*L + phase;
                y[yi] = apply_filter(&x[j], phase);
            }
        }

        // push end of buffer
        const int M1 = _max(N-K, M0);
        push_values(&x[M1], N-M1);
    }

private:
    inline
    void push_value(T x) {
        for (int i = 0; i < (K-1); i++) {
            xn[i] = xn[i+1];
        }
        xn[K-1] = x;
    }

    inline
    void push_values(const T* x, const int N) {
        const int M = K-N;
        for (int i = 0; i < M; i++) {
            xn[i] = xn[i+N];
        }
        for (int i = M, j = 0; i < K; i++, j++) {
            xn[i] = x[j];
        }
    }

    inline
    T apply_filter(const T* x, const int phase) {
        auto* b0 = &b[phase*K];
        T y; 
        y = 0;
        for (int i = 0; i < K; i++) {
            y += x[i] * b0[i];
        }
        return y;
    }
};

#undef _min
#undef _max

#include "simd/f32_cum_mul.h"
template <>
inline
float PolyphaseDownsampler<float>::apply_filter(const float* x) {
    return f32_cum_mul_auto(x, b.data(), NN);
}

#include "simd/c32_f32_cum_mul.h"
template <>
inline
std::complex<float> PolyphaseDownsampler<std::complex<float>>::apply_filter(const std::complex<float>* x) {
    return c32_f32_cum_mul_auto(x, b.data(), NN);
}