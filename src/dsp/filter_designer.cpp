#define _USE_MATH_DEFINES
#include <cmath>
#include <assert.h>
#include <complex>

#include "filter_designer.h"

constexpr 
float PI = (float)M_PI;

inline
auto get_phasor(float x) {
    return std::complex<float>(
        std::cos(x),
        std::sin(x)
    );
}

inline
float sinc(float x) {
    if (std::abs(x) <= 1e-6f) {
        return 1.0f;
    }
    return std::sin(PI*x)/(PI*x);
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

// Apply pre-warping to compensate for non-linear frequency mapping of bilinear transform
float prewarp_normalised_frequency(const float Kd) {
    // Consider the bilinear transform
    // z = exp(sT)
    // s = 1/T ln(z)
    // s ≈ 2/T * (z-1)/(z+1)

    // Wd = actual discrete filter frequency, Wa = the analog frequency we use in raw bilinear transform

    // Hd(z) = Ha(s)
    // Hd(exp(jwdT)) = Ha(2/T * (exp(jwdT)-1)/(exp(jwdT)+1))
    // exp(jwdT) = 2/T * (exp(jwdT)-1)/(exp(jwdT)+1)
    //           = 2/T * exp(jwdT/2)/exp(jwdT/2) * [exp(jwdT/2) - exp(-jwdT/2)]/[exp(jwdT/2) + exp(-jwdT/2)]
    //           = 2/T * j * sin(wdT/2)cos(wdT/2)
    //           = 2/T * j * tan(wd*T/2)
    // wa = 2/T * tan(wd*T/2)
    
    // Wa = 2/T tan(Wd T/2)
    // 2*pi*Fa = 2*Fs tan[(2*pi*Fd)/(2*Fs)]
    // Fa/(Fs/2) = 2/pi tan[pi/2 * Fd/(Fs/2)]
    // Ka = 2/pi tan(pi/2 * Kd)
    const float Ka = 2.0f/PI * std::tan(PI/2.0f * Kd);
    return Ka;
}

// For FIR filter design we consider two things
// 1. The ideal brickwall response given by an infinite sinc(t) function in time domain
// 2. The response due to the windowing function
//
// Ideal brickwall response is given as
// sinc(t) <--> Fourier <--> H_brickwall(f)
//
// The rectangular window is used whenever we limit the number of samples of the ideal brickwall in time domain
// h_rect(t) <--> Fourier <--> sinc(f)
//
// The resulting filter is a timewise multiplication of the ideal brickwall and window function
// Using Fourier transform properties, this is convolution in the frequency domain
// sinc(t)*h_rect(t) <--> Fourier <--> H_brickwall(f) <convolve> sinc(f)
// This produces sidelobes of high magnitude
//
// Thus we pick a better window function apart from a rectangular window
// The Hamming window in frequency domain has sidelobes of lower magnitude

void create_fir_lpf(float* b, const int N, const float k, const window_func_t window) {
    assert(b != NULL);
    assert(N > 0);
    assert(k < 1.0f);
    assert(k > 0.0f);

    auto _b = ReverseArray(b, N);
    const float M = (float)(N-1);
    for (int i = 0; i < N; i++) {
        const float t0 = 2.0f*PI*(float)(i)/M;
        const float t1 = (float)i - M/2.0f;

        const float h_window = window(t0);
        // Brickwall response in frequency can be scaled
        // Using the scaling property of the Fourier transform
        // h(at) <--> Fourier <--> H(w/a)/|a| 
        // |a|*h(at) <--> Fourier <--> H(w/a)
        // Increasing a makes H(w/a) wider
        // sinc(t) <--> Fourier <--> rect(w)
        // k*sinc(k*t) <--> Fourier <--> rect(w/k)
        const float h_filter = k*sinc(k*t1);
        _b[i] = h_window*h_filter; 
    }
}

void create_fir_hpf(float* b, const int N, const float k, const window_func_t window) {
    assert(b != NULL);
    assert(N > 0);
    assert(k < 1.0f);
    assert(k > 0.0f);

    auto _b = ReverseArray(b, N);
    const float M = (float)(N-1);
    for (int i = 0; i < N; i++) {
        const float t0 = 2*PI*(float)(i)/M;
        const float t1 = (float)(i) - M/2.0f;

        const float h_window = window(t0);
        // High pass filter is H_high(f) = 1-H_low(f) 
        // H_high(f) <--> Fourier <--> h_high(t)
        // H_high(f) <--> Fourier <--> sinc(t) - h_low(t)
        // h_high(t) = sinc(t) - h_low(t)
        const float h_filter = sinc(t1) - k*sinc(k*t1);
        _b[i] = h_window*h_filter; 
    }
}

void create_fir_bpf(float* b, const int N, const float k1, const float k2, const window_func_t window) {
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

        const float h_window = window(t0);
        // Bandpass filter is the subtraction of two low pass filters
        // H_band(f) = H_low_upper(f) - H_low_lower(f)
        // H_band(f) <--> Fourier <--> h_band(t)
        // H_band(f) <--> Fourier <--> h_low_upper(t) - h_low_lower(t)
        // h_band(t) = h_low_upper(t) - h_low_lower(t)
        const float h_filter = k2*sinc(k2*t1) - k1*sinc(k1*t1);
        _b[i] = h_window*h_filter; 
    }
}


void create_iir_single_pole_lpf(float* b, float* a, const float k) {
    assert(b != NULL);
    assert(a != NULL);
    assert(k < 1.0f);
    assert(k > 0.0f);

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

    const float k_warp = prewarp_normalised_frequency(k);
    const float A = 1.0f/(PI*k_warp);
    const float B0 = 1.0f + 2.0f*A;
    const float B1 = 1.0f - 2.0f*A;
    const float b0 = 1.0f/B0;
    const float a0 = B1/B0;

    const int N = TOTAL_TAPS_IIR_SINGLE_POLE_LPF;
    auto _b = ReverseArray(b, N);
    auto _a = ReverseArray(a, N);

    _b[0] = b0;
    _b[1] = b0;
    _a[0] = 1.0f;
    _a[1] = -a0;
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

    // Use exact z transform
    // z = exp(sT), T = 1/Fs
    // z = exp(j*w0*Ts)
    // Let wn = w0*Ts = 2*pi*Fc/Fs
    // k = Fc/(Fs/2) = 2Fc/Fs 
    // wn = pi*k

    // Use pole placement technique in z-domain
    // Place zeros on unit circle at ±wn
    // Place poles slightly inside unit circle at ±wn
    // The close the poles are to the unit circle, the tighter the notch response
    // HINT: Visualise the length of the vector [z-exp(jwn)] as z moves along the unit circle
    // H(z) = [(z-exp(jwn))(z-exp(-jwn))] / [(z-r*exp(jwn))(z-r*exp(-jwn))]
    // H(z) = [z^2 - 2cos(wn)*z + 1] / [z^2 - 2rcos(wn)*z + r^2] 
    // H(z) = [1 - 2cos(wn)*z^-1 + z^-2] / [1 - 2rcos(wn)*z^-1 + (r^2)*z^-2]

    const float wn = PI*k;
    const float a0 = 2.0f*std::cos(wn);
    const float r2 = r*r;

    static auto H_z = [k,r](float k_z) {
        // z = exp(sT) = exp(j*pi*k), k = Fc/(Fs/2)
        const auto z = get_phasor(PI*k_z);
        const auto z0 = get_phasor(+PI*k);
        const auto z1 = get_phasor(-PI*k);
        return ((z-z0)*(z-z1))/((z-r*z0)*(z-r*z1));
    };

    // Find a frequency which we normalise the magnitude to
    const float k_z = (k > 0.5f) ? 0.0f : 1.0f;
    const float K = 1.0f/std::abs(H_z(k_z));

    _b[0] = 1.0f;
    _b[1] = -a0;
    _b[2] = 1.0f;

    _a[0] = 1.0f;
    _a[1] = a0*r;
    _a[2] = -r2;

    for (int i = 0; i < N; i++) {
        _b[i] = K*_b[i];
    }
}

void create_iir_peak_1_filter(float* b, float* a, const float k, const float r) {
    assert(b != NULL);
    assert(a != NULL);
    assert(k < 1.0f);
    assert(k > 0.0f);
    assert(r < 1.0f);
    assert(r > 0.0f);

    const int N = TOTAL_TAPS_IIR_SECOND_ORDER_PEAK_FILTER;
    auto _b = ReverseArray(b, N);
    auto _a = ReverseArray(a, N);

    // Use exact z transform
    // z = exp(sT), T = 1/Fs
    // z = exp(j*w0*Ts)
    // Let wn = w0*Ts = 2*pi*Fc/Fs
    // k = Fc/(Fs/2) = 2Fc/Fs 
    // wn = pi*k

    // Use pole placement
    // H(z) = 1/[(z-r*exp(jwn))*(z-r*exp(-jwn))]
    // H(z) = 1/[z^2 - 2*r*cos(wn)*z + (r^2)]
    // H(z) = (z^-2)/[1 - 2*r*cos(wn)*z^-1 + (r^2)*z^-2]

    const float wn = PI*k;
    const float a0 = 2.0f*std::cos(wn);
    const float r2 = r*r;

    static auto H_z = [k,r](float k_z) {
        // z = exp(sT) = exp(j*pi*k), k = Fc/(Fs/2)
        const auto z = get_phasor(PI*k_z);
        const auto z0 = get_phasor(+PI*k);
        const auto z1 = get_phasor(-PI*k);
        return 1.0f/((z-r*z0)*(z-r*z1));
    };

    // Normalise to the peak frequency
    const float K = 1.0f/std::abs(H_z(k));

    _b[0] = 0.0f;
    _b[1] = 0.0f;
    _b[2] = 1.0f;

    _a[0] = 1.0f;
    _a[1] = r*a0;
    _a[2] = -r2;

    for (int i = 0; i < N; i++) {
        _b[i] = K*_b[i];
    }
}

void create_iir_peak_2_filter(float* b, float* a, const float k, const float r, const float A_db) {
    assert(b != NULL);
    assert(a != NULL);
    assert(k < 1.0f);
    assert(k > 0.0f);
    assert(r < 1.0f);
    assert(r > 0.0f);

    const int N = TOTAL_TAPS_IIR_SECOND_ORDER_PEAK_FILTER;
    auto _b = ReverseArray(b, N);
    auto _a = ReverseArray(a, N);

    // Use exact z transform
    // z = exp(sT), T = 1/Fs
    // z = exp(j*w0*Ts)
    // Let wn = w0*Ts = 2*pi*Fc/Fs
    // k = Fc/(Fs/2) = 2Fc/Fs 
    // wn = pi*k

    // Use zero and pole placement
    // H(z) = [(z-r0*exp(jwn))*(z-r0*exp(-jwn))]/[(z-r1*exp(jwn))*(z-r1*exp(-jwn))]
    // H(z) = [z^2 - 2*r0*cos(wn)*z + (r0^2)]/[z^2 - 2*r1*cos(wn)*z + (r1^2)]
    // H(z) = [1 - 2*r0*cos(wn)*z^-1 + (r0^2)*z^-2]/[1 - 2*r1*cos(wn)*z^-1 + (r1^2)*z^-2]

    const float A = std::pow(10.0f, A_db/20.0f);
    const float rc = 1.0f-r;
    const float rc_scale = rc*2.0f;
    const float r0 = 1.0f - rc_scale; 
    const float r1 = 1.0f - rc_scale/A;

    const float wn = PI*k;
    const float a0 = 2.0f*std::cos(wn);

    static auto H_z = [k,r0,r1](float k_z) {
        // z = exp(sT) = exp(j*pi*k), k = Fc/(Fs/2)
        const auto z = get_phasor(PI*k_z);
        const auto z0 = get_phasor(+PI*k);
        const auto z1 = get_phasor(-PI*k);
        return (z-r0*z0)*(z-r0*z1)/((z-r1*z0)*(z-r1*z1));
    };

    // Normalise to the peak frequency
    const float K = 1.0f/std::abs(H_z(k));

    _b[0] = 1.0f;
    _b[1] = -r0*a0;
    _b[2] = r0*r0;

    _a[0] = 1.0f;
    _a[1] = r1*a0;
    _a[2] = -r1*r1;

    for (int i = 0; i < N; i++) {
        _b[i] = K*_b[i];
    }
}

void create_fir_hilbert(float* b, const int N) {
    assert(b != NULL);
    assert(N > 0);

    auto _b = ReverseArray(b, N);

    // Given by the non-causal equation 
    // n != 0 : 2/(n*pi) * sin(n*pi/2)^2 
    // n  = 0 : 0
    const int M = (N-1)/2;
    for (int i = 0; i < N; i++) {
        const int n = i-M;
        const bool is_even = (n % 2) == 0;
        _b[i] = is_even ? 0.0f : 2.0f/(PI*(float)n);
    }
}