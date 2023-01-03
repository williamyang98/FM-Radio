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
    assert(k < 1.0f);
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
    assert(k < 1.0f);
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
    assert(k2 > k1);
    assert(k2 < 1.0f);
    assert(k1 < 1.0f);
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
    // Creates a single pole iir filter using bilinear transform
    // Apply pre-warping to compensate for non-linear frequency mapping of bilinear transform
    // Wd = actual discrete filter frequency, Wa = the analog frequency we use in raw bilinear transform
    // Wa = 2/T tan(Wd T/2)
    // 2*pi*Fa = 2*Fs tan[(2*pi*Fd)/(2*Fs)]
    // Fa/(Fs/2) = 2/pi tan[pi/2 * Fd/(Fs/2)]
    // Ka = 2/pi tan(pi/2 * Kd)
    const float k_warp = 2.0f/PI * std::tan(PI/2.0f * k);

    // k' = Tc/Ts
    // Tc = 1/(2*pi*fc)
    // k' = Fs/fc * 1/(2*pi)
    // k' = (Fs/2)/fc * 1/pi
    // k' = 1/(k*pi)
    const float z = 1.0f/(k_warp*PI);

    const float b0 = 1.0f/(1.0f+2.0f*z);
    const float a0 = (1.0f-2.0f*z)/(1.0f+2.0f*z);

    const int N = 2;
    auto _b = ReverseArray(b, N);
    auto _a = ReverseArray(a, N);

    _b[0] = b0;
    _b[1] = b0;
    _a[0] = 1.0f;
    _a[1] = a0;
}

void create_iir_notch_filter(float* b, float* a, const float k, const float r) {
    const float w0 = PI*k;
    const float a0 = 2.0f*std::cos(w0);

    const int N = 3;
    auto _b = ReverseArray(b, N);
    auto _a = ReverseArray(a, N);

    _b[0] = 1.0f;
    _b[1] = -a0;
    _b[2] = 1.0f;

    _a[0] = 1.0f;
    _a[1] = -a0*r;
    _a[2] = r*r;
}

void create_iir_peak_filter(float* b, float* a, const float k, const float r) {
    const float w0 = PI*k;
    const float a0 = 2.0f*std::cos(w0);

    const int N = 3;
    auto _b = ReverseArray(b, N);
    auto _a = ReverseArray(a, N);

    _b[0] = 0.0f;
    _b[1] = a0*(1.0f-r);
    _b[2] = r*r - 1.0f;

    _a[0] = 1.0f;
    _a[1] = -a0*r;
    _a[2] = r*r;
}

void create_fir_hilbert(float* b, const int N) {
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