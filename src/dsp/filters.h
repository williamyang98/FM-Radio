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
    FIR_Filter(const float* _b, const int _K) 
    : K(_K), b(_K), xn(_K)
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
        const int M0 = _min(K-1, N);

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
        const int M1 = _max(N-K, M0);
        push_values(&x[M1], N-M1);
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

template <typename T>
class IIR_Filter 
{
public:
    const int K;
    int Ki;
    float* b; 
    float* a; 
    T* xn;
    T* yn;
public:
    IIR_Filter(const float* _b, const float* _a, const int _K) 
    : K(_K)
    {
        b = new float[K];
        a = new float[K];
        for (int i = 0; i < K; i++) {
            b[i] = _b[i];
            a[i] = _a[i];
        }
        xn = new T[K]{};
        yn = new T[K]{};
        Ki = 0;
    }

    void process(const T* x, T* y, const int N) {
        for (int i = 0; i < N; i++) {
            xn[Ki] = x[i];
            y[i] = xn[Ki]*b[0];
            for (int j = 1; j < K; j++) {
                const int k0 = (K+Ki-j)%K;
                y[i] += (xn[k0]*b[j] - yn[k0]*a[j]);
            }
            yn[Ki] = y[i];
            Ki = (Ki+1)%K;
        }
    }

    ~IIR_Filter() {
        delete [] b;
        delete [] a;
        delete [] xn;
        delete [] yn;
    }
};

template <typename T>
class AGC_Filter 
{
public:
    float target_power = 1.0f;
    float current_gain = 0.1f;
    float beta = 0.2f;
    void process(const T* x, T* y, const int N) {
        const float avg_power = calculate_average_power(x, N);
        const float target_gain = std::sqrt(target_power/avg_power);
        current_gain = current_gain + beta*(target_gain - current_gain);
        for (int i = 0; i < N; i++) {
            y[i] = current_gain*x[i];
        }
    }
private:
    float calculate_average_power(const T* x, const int N) {
        float avg_power = 0.0f;
        for (int i = 0; i < N; i++) {
            const float I = x[i].real();
            const float Q = x[i].imag();
            avg_power += (I*I + Q*Q);
        }
        avg_power /= (float)N;
        return avg_power;
    }
};

template <typename T>
class Integrator_Block 
{
public:
    float KTs = 1.0f;
    T yn = 0.0f;
    T process(T x) {
        T y = KTs*x + yn; 
        yn = y;
        return y;
    }
};

template <typename T>
class Delay_Line
{
public:
    T* xn;
    const int K;
    int Ki;
public:
    Delay_Line(const int _K, T v = 0)
    : K(_K) 
    {
        xn = new T[K]{v};
        Ki = 0;
    }
    ~Delay_Line() {
        delete [] xn;
    }
    void process(const T* x, T* y, const int N) {
        for (int i = 0; i < N; i++) {
            xn[Ki] = x[i];
            y[i] = xn[(Ki+K-1)%K];
            Ki = (Ki+1)%K;
        }
    }
};

#undef _min
#undef _max