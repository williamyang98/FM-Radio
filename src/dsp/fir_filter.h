#pragma once
#include "utility/aligned_vector.h"

#define _min(A,B) (A > B) ? B : A
#define _max(A,B) (A > B) ? A : B

template <typename T>
class FIR_Filter 
{
public:
    const int K;
    AlignedVector<float> b;
    AlignedVector<T> xn;
    AlignedVector<T> tmp;
public:
    FIR_Filter(const float* _b, const int _K) 
    : K(_K), b(_K), xn(_K), tmp(_K)
    {
        // time reverse the filter
        for (int i = 0; i < K; i++) {
            b[i] = _b[(K-1)-i];
        }

        for (int i = 0; i < K; i++) {
            xn[i] = 0;
        }
    }

    void process(const T* x, T* y, const int N) {
        const int M0 = _min(K-1, N);  // head
        const int M1 = _max(N-K, M0); // tail

        // save the tail which is overridden
        for (int i = M1, j = 0; i < N; i++, j++) {
            tmp[j] = x[i];
        }

        // inplace math 
        // NOTE: We do this first and in reverse so that we don't modify values
        //       if x and y are the same buffer
        // Forward loop: for (int i = M0, j = 0; i < N; i++, j++)
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
        T y = 0;
        for (int i = 0; i < K; i++) {
            y += (x[i] * b[i]);
        }
        return y;
    }
};

#undef _min
#undef _max