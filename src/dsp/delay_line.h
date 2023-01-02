#pragma once
#include "utility/aligned_vector.h"

#define _min(A,B) (A > B) ? B : A
#define _max(A,B) (A > B) ? A : B

template <typename T>
class Delay_Line
{
public:
    AlignedVector<T> xn;
    AlignedVector<T> tmp;
    const int K;
public:
    Delay_Line(const int _K, T v = 0)
    : K(_K), 
      xn(K),
      tmp(K)
    {
        for (int i = 0; i < K; i++) {
            xn[i] = v;
        }
    }
    void process(const T* x, T* y, const int N) {
        const int M0 = _min(K, N);  // head
        const int M1 = _max(N-K, 0); // tail

        // save the tail which is overridden
        for (int i = M1, j = 0; i < N; i++, j++) {
            tmp[j] = x[i];
        }

        // inplace shift
        // NOTE: We do this first and in reverse so that we don't modify values
        //       if x and y are the same buffer
        // Forward loop: for (int i = M0, j = 0; i < N; i++, j++)
        for (int i = (N-1), j = (N-1)-K; i >= M0; i--, j--) {
            y[i]  = x[j];
        }

        // continue from previous block
        for (int i = 0; i < M0; i++) {
            y[i] = xn[i];
        }

        // push end of buffer
        push_values(tmp.data(), N-M1);
    }

private:
    void push_values(const T* x, const int N) {
        const int M = K-N;
        for (int i = 0; i < M; i++) {
            xn[i] = xn[i+N];
        }
        for (int i = M, j = 0; i < K; i++, j++) {
            xn[i] = x[j];
        }
    }
};

#undef _min
#undef _max