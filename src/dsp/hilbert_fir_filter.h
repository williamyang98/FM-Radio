#pragma once

#include "fir_filter.h"
#include "filter_designer.h"
#include <complex>

#define _min(A,B) (A > B) ? B : A
#define _max(A,B) (A > B) ? A : B

// Implement the Hilbert transform using an FIR filter
template <typename T>
class Hilbert_FIR_Filter: public FIR_Filter<T>
{
public:
    Hilbert_FIR_Filter(const int _K)
    : FIR_Filter(_K)
    {
        create_fir_hilbert(b.data(), K);
    }

    // Generate the quadrature component
    // Real component is delayed by floor(K/2) samples to match delay of hilbert filter
    // Imag component contains the Hilbert transform delayed by floor(K/2) samples
    void process(const T* x, std::complex<T>* y, const int N) {
        const int M = (K-1)/2;
        const int M0 = _min(K-1, N);
        const int M1 = _max(N-K, M0);

        // continue from previous block
        for (int i = 0; i < M0; i++) {
            push_value(x[i]);
            const T imag = apply_filter(xn.data());
            y[i] = { xn[M], imag };
        }

        // push end of buffer
        push_values(&x[M1], N-M1);

        // inplace math
        for (int i = M0, j = 0; i < N; i++, j++) {
            const T imag = apply_filter(&x[j]);
            y[i] = { x[j+M], imag };
        }
    }
};

#undef _min
#undef _max