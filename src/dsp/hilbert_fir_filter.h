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
private:
    using Base = FIR_Filter<T>;
public:
    Hilbert_FIR_Filter(const int _K)
    : FIR_Filter<T>(_K)
    {
        create_fir_hilbert(Base::b.data(), Base::K);
    }

    // Generate the quadrature component
    // Real component is delayed by floor(K/2) samples to match delay of hilbert filter
    // Imag component contains the Hilbert transform delayed by floor(K/2) samples
    void process(const T* x, std::complex<T>* y, const int N) {
        const int M = (Base::K-1)/2;
        const int M0 = _min(Base::K-1, N);
        const int M1 = _max(N-Base::K, M0);

        // continue from previous block
        for (int i = 0; i < M0; i++) {
            Base::push_value(x[i]);
            const T imag = Base::apply_filter(Base::xn.data());
            y[i] = { Base::xn[M], imag };
        }

        // push end of buffer
        Base::push_values(&x[M1], N-M1);

        // inplace math
        for (int i = M0, j = 0; i < N; i++, j++) {
            const T imag = Base::apply_filter(&x[j]);
            y[i] = { x[j+M], imag };
        }
    }
};

#undef _min
#undef _max