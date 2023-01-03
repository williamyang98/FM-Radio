#define _USE_MATH_DEFINES
#include <cmath>
#include <assert.h>

#include "filter_designer.h"

constexpr 
float PI = (float)M_PI;

inline
float sinc(float x) {
    if (std::abs(x) <= 1e-6f) {
        return 1.0f;
    }
    return std::sin(PI*x)/(PI*x);
}

inline
float calc_hamming_window(float x) {
    return 0.54f - 0.46f*std::cos(x);
}

template <typename T>
struct ReverseArray 
{
private:
    T* x;
    const int N;
public:
    ReverseArray(T* _x, const int _N)
    : x(_x), N(_N) {}
    T& operator[](const int i) {
        return x[(N-1)-i];
    }
};

void create_fir_lpf(float* b, const int N, const float k) {
    assert(b != NULL);
    assert(N > 0);
    assert(k < 1.0f);
    assert(k > 0.0f);

    auto _b = ReverseArray(b, N);
    const float M = (float)(N-1);
    for (int i = 0; i < N; i++) {
        const float t0 = 2.0f*PI*(float)(i)/M;
        const float t1 = (float)i - M/2.0f;

        const float h_window = calc_hamming_window(t0);
        const float h_filter = k*sinc(k*t1);
        _b[i] = h_window*h_filter; 
    }
}

void create_fir_hpf(float* b, const int N, const float k) {
    assert(b != NULL);
    assert(N > 0);
    assert(k < 1.0f);
    assert(k > 0.0f);

    auto _b = ReverseArray(b, N);
    const float M = (float)(N-1);
    for (int i = 0; i < N; i++) {
        const float t0 = 2*PI*(float)(i)/M;
        const float t1 = (float)(i) - M/2.0f;

        const float h_window = calc_hamming_window(t0);
        const float h_filter = sinc(t1) - k*sinc(k*t1);
        _b[i] = h_window*h_filter; 
    }
}

void create_fir_bpf(float* b, const int N, const float k1, const float k2) {
    assert(b != NULL);
    assert(N > 0);
    assert(k2 > k1);
    assert(k2 < 1.0f);
    assert(k1 < 1.0f);
    assert(k1 > 0.0f);
    assert(k2 > 0.0f);

    auto _b = ReverseArray(b, N);
    const float M = (float)N-1;
    for (int i = 0; i < N; i++) {
        const float t0 = 2*PI*(float)(i)/M;
        const float t1 = (float)(i) - M/2.0f;

        const float h_window = calc_hamming_window(t0);
        const float h_filter = k2*sinc(k2*t1) - k1*sinc(k1*t1);
        _b[i] = h_window*h_filter; 
    }
}


void create_iir_single_pole_lpf(float* b, float* a, const float k) {
    assert(b != NULL);
    assert(a != NULL);
    assert(k < 1.0f);
    assert(k > 0.0f);

    // Apply pre-warping to compensate for non-linear frequency mapping of bilinear transform
    // Wd = actual discrete filter frequency, Wa = the analog frequency we use in raw bilinear transform
    // Wa = 2/T tan(Wd T/2)
    // 2*pi*Fa = 2*Fs tan[(2*pi*Fd)/(2*Fs)]
    // Fa/(Fs/2) = 2/pi tan[pi/2 * Fd/(Fs/2)]
    // Ka = 2/pi tan(pi/2 * Kd)
    const float k_warp = 2.0f/PI * std::tan(PI/2.0f * k);

    // Transfer function of analog single pole low pass filter
    // H(s) = 1/(1 + Tc*s)
    // Tc = 1/(2*pi*fc) = 1/wc
    
    // Using bilinear transform
    // s = 2/Ts * (1-z^-1) / (1+z^-1)
    // Y(s) = X(s) * 1/(1+Tc*s)
    // Y(s) * (1 + Tc*s) = X(s)
    // Y(z) * (1 + 2*Tc/Ts * (1-z^-1)/(1+z^-1) ) = X(z)

    // Let A = Tc/Ts
    // A = Fs/(2*pi*fc) = Fs/(2Fc) * 1/pi
    // k = Fc/(Fs/2) = 2Fc/Fs
    // A = 1/(pi*k)

    // Y(z) * [(1+z^-1) + 2A*(1-z^-1)] = X(z) * (1+z^-1)
    // Y(z) * [(1+2A) + (1-2A)*z^-1] = X(z) * (1+z^-1)
    // (1+2A)*y[n] + (1-2A)*y[n-1] = x[n] + x[n-1]
    // y[n] = 1/(1+2A)*x[n] +      1/(1+2A)*x[n-1]
    //                      - (1-2A)/(1+2A)*y[n-1]

    const float A = 1.0f/(PI*k_warp);
    const float B0 = 1.0f + 2.0f*A;
    const float B1 = 1.0f - 2.0f*A;
    const float b0 = 1.0f/B0;
    const float a0 = -B1/B0;

    const int N = TOTAL_TAPS_IIR_SINGLE_POLE_LPF;
    auto _b = ReverseArray(b, N);
    auto _a = ReverseArray(a, N);

    _b[0] = b0;
    _b[1] = b0;
    _a[0] = 1.0f;
    _a[1] = a0;
}

void create_iir_notch_filter(float* b, float* a, const float k, const float r) {
    assert(b != NULL);
    assert(a != NULL);
    assert(k < 1.0f);
    assert(k > 0.0f);
    assert(r < 1.0f);
    assert(r > 0.0f);


    const int N = TOTAL_TAPS_IIR_SECOND_ORDER_NOTCH_FILTER;
    auto _b = ReverseArray(b, N);
    auto _a = ReverseArray(a, N);

    // https://dsp.stackexchange.com/a/31030
    // wn = 2*pi*Fc/Fs
    // k = Fc/(Fs/2) = 2*Fc/Fs
    // wn = pi*k
    const float wn = PI*k;
    const float a0 = 2.0f*std::cos(wn);
    const float K = (1.0f + r*a0 + r*r) / (2.0f + a0);

    _b[0] = K;
    _b[1] = -K*a0;
    _b[2] = K;

    _a[0] = 1.0f;
    _a[1] = a0*r;
    _a[2] = -r*r;
}

void create_iir_peak_filter(float* b, float* a, const float k, const float r) {
    assert(b != NULL);
    assert(a != NULL);
    assert(k < 1.0f);
    assert(k > 0.0f);
    assert(r < 1.0f);
    assert(r > 0.0f);

    const int N = TOTAL_TAPS_IIR_SECOND_ORDER_PEAK_FILTER;
    auto _b = ReverseArray(b, N);
    auto _a = ReverseArray(a, N);

    // wn = 2*pi*Fc/Fs
    // k = Fc/(Fs/2) = 2*Fc/Fs
    // wn = pi*k
    const float wn = PI*k;
    const float a0 = 2.0f*std::cos(wn);

    _b[0] = 0.0f;
    _b[1] = a0*(1.0f-r);
    _b[2] = r*r - 1.0f;

    _a[0] = 1.0f;
    _a[1] = a0*r;
    _a[2] = -r*r;
}

void create_fir_hilbert(float* b, const int N) {
    assert(b != NULL);
    assert(N > 0);

    // Given by the non-causal equation 
    // n != 0 : 2/(n*pi) * sin(n*pi/2)^2 
    // n  = 0 : 0
    const int M = (N-1)/2;
    for (int i = 0; i < N; i++) {
        const int n = i-M;
        if (n == 0) {
            b[i] = 0.0f;
            continue;
        }

        // NOTE: time reverse coefficients here
        const float _n = -(float)n;  
        const float a0 = std::sin(_n*PI/2.0f);
        b[i] = 2.0f/(_n*PI) * a0 * a0;
    }
}