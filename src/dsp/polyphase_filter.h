#pragma once

template <typename T>
class PolyphaseUpsampler 
{
private:
    const int L;
    const int K;
    const int NN;
    float* b;
    T* xn;
    int Ki;
public:
    // b = FIR filter coefficients of length L*K
    // L = upsampling factor and total phases
    // K = number of coefficients per phase
    PolyphaseUpsampler(const float *_b, const int _L, const int _K) 
    : L(_L), K(_K), NN(_L*_K) {
        b = new float[NN];
        for (int i = 0 ; i < NN; i++) {
            b[i] = _b[i] * (float)L;
        }

        xn = new T[K]{0};
        Ki = 0;
    }

    ~PolyphaseUpsampler() {
        delete [] b;
        delete [] xn;
    }

    // E.g. L = 4, K = 2, L*K = 8
    //  0  0  0 x0  0  0  0 x1  0  0  0 x2 
    // b7 b6 b5 b4 b3 b2 b1 b0          => y0
    //    b7 b6 b5 b4 b3 b2 b1 b0       => y1
    //       b7 b6 b5 b4 b3 b2 b1 b0    => y2
    //          b7 b6 b5 b4 b3 b2 b1 b0 => y3
    // N = process N input samples
    void process(const T* x, T* y, const int N) {
        for (int i = 0; i < N; i++) {
            xn[Ki] = x[i];

            for (int j = 0; j < L; j++) {
                const int yi = i*L + j;
                y[yi] = 0;
                for (int k = 0; k < K; k++) {
                    const int xi = (K+Ki-k)%K;
                    y[yi] += xn[xi] * b[j + k*L];
                }
            }
            Ki = (Ki+1)%K;
        }
    }
};

template <typename T>
class PolyphaseDownsampler 
{
private:
    const int M;
    const int K;
    const int NN;
    float* b;
    T* xn;
    int Ki;
public:
    // b = FIR filter with M*K coefficients
    // M = downsampling factor and total phases 
    // K = total coefficients per phase
    PolyphaseDownsampler(const float *_b, const int _M, const int _K) 
    : M(_M), K(_K), NN(_M*_K) {
        b = new float[NN];
        for (int i = 0 ; i < NN; i++) {
            b[i] = _b[i];
        }

        xn = new T[NN]{0};
        Ki = 0;
    }

    ~PolyphaseDownsampler() {
        delete [] b;
        delete [] xn;
    }

    // E.g. M = 3, K = 2, M*K = 6
    // x8 x7 x6 x5 x4 x3 x2 x1 x0
    // b5 b4 b3 b2 b1 b0          => y0
    //    b5 b4 b3 b2 b1 b0       => ... (downsampling ignores these)
    //       b5 b4 b3 b2 b1 b0    => ...
    //          b5 b4 b3 b2 b1 b0 => y1
    // N = produce N output samples
    void process(const T* x, T* y, const int N) {
        for (int xi = 0, yi = 0; yi < N; yi++, xi+=M) {
            // get M samples and distribute along phases
            for (int i = 0; i < M; i++) {
                Ki = (Ki+1)%NN;
                xn[Ki] = x[xi+i];
            }

            y[yi] = 0;
            for (int i = 0; i < NN; i++) {
                const int xni = (NN+Ki-i) % NN;
                y[yi] += xn[xni] * b[i];
            }
        }
    }
};


