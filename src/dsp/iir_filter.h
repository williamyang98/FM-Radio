#pragma once
#include "utility/aligned_vector.h"

template <typename T>
class IIR_Filter
{
private:
    const int K;
    AlignedVector<float> b;
    AlignedVector<float> a;
    AlignedVector<T> xn;
    AlignedVector<T> yn;
public:
    float* get_b() const { return b.data(); }
    float* get_a() const { return a.data(); }
    int    get_K() const { return K; }
public:
    IIR_Filter(const int _K) 
    : K(_K),
      b(K), a(K),
      xn(K), yn(K)
    {
        for (int i = 0; i < K; i++) {
            b[i] = 0;
            a[i] = 0; 
            xn[i] = 0;
            yn[i] = 0;
        }

        // Direct form I means that this should just be 1.0f
        // H(z) = (b0 + b1*z^-1 + b2*z^-2 + ...) / (1 - a1*z^-1 - a2*z^-2 - ...) 
        // H(z) = Y(z)/X(z)
        // y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] + ...
        //                + a1*y[n-1] + a2*y[n-2] + ...
        // Therefore we ignore the calculation of a0*y[n] which is not used
        yn[K-1] = 0;
    }

    void process(const T* x, T* y, const int N) {
        for (int i = 0; i < N; i++) {
            push_x(x[i]);
            y[i] = apply_filter();
            push_y(y[i]);
        }
    }
private:
    void push_x(T x) {
        for (int i = 0; i < (K-1); i++) {
            xn[i] = xn[i+1];
        }
        xn[K-1] = x;
    }

    void push_y(T y) {
        for (int i = 0; i < (K-2); i++) {
            yn[i] = yn[i+1];
        }
        yn[K-2] = y;
    }

    T apply_filter() {
        T y; 
        y = 0;
        for (int i = 0; i < K; i++) {
            y += (xn[i]*b[i] + yn[i]*a[i]);
        }
        return y;
    }
};