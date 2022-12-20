#define _USE_MATH_DEFINES
#include <cmath>
#include <assert.h>

#include "filter_designer.h"

const float PI = (float)M_PI;

float sinc(float x) {
    if (std::abs(x) <= 1e-6f) {
        return 1.0f;
    }
    return std::sin(PI*x)/(PI*x);
}

float calc_hamming_window(float x) {
    return 0.54f - 0.46f*std::cos(x);
}

// k = Fc/(Fs/2)
FIR_Filter_Res* create_fir_lpf(const float k, const int M) {
    assert(k < 1.0f);
    auto f = new FIR_Filter_Res(M+1);
    for (int n = 0; n <= M; n++) {
        const float h_window = calc_hamming_window(2*PI*(float)(n)/(float)(M));
        const float n_shift = (float)(n) - (float)(M)/2.0f;
        const float h_filter = k*sinc(k*n_shift);
        f->b[n] = h_window*h_filter; 
    }
    return f;
}

// k = Fc/(Fs/2)
FIR_Filter_Res* create_fir_hpf(const float k, const int M) {
    assert(k < 1.0f);
    auto f = new FIR_Filter_Res(M+1);
    for (int n = 0; n <= M; n++) {
        const float h_window = calc_hamming_window(2*PI*(float)(n)/(float)(M));
        const float n_shift = (float)(n) - (float)(M)/2.0f;
        const float h_filter = sinc(n_shift) - k*sinc(k*n_shift);
        f->b[n] = h_window*h_filter; 
    }
    return f;
}

// k = Fc/(Fs/2)
FIR_Filter_Res* create_fir_bpf(const float k1, const float k2, const int M) {
    assert(k2 > k1);
    assert(k2 < 1.0f);
    assert(k1 < 1.0f);

    auto f = new FIR_Filter_Res(M+1);
    for (int n = 0; n <= M; n++) {
        const float h_window = calc_hamming_window(2*PI*(float)(n)/(float)(M));
        const float n_shift = (float)(n) - (float)(M)/2.0f;
        const float h_filter = k2*sinc(k2*n_shift) - k1*sinc(k1*n_shift);
        f->b[n] = h_window*h_filter; 
    }
    return f;
}


// Creates a single pole iir filter using bilinear transform
// k = Fc/(Fs/2)
IIR_Filter_Res* create_iir_single_pole_lpf(const float k) {
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

    auto f = new IIR_Filter_Res(2);

    const float b = 1.0f/(1.0f+2.0f*z);
    const float a = (1.0f-2.0f*z)/(1.0f+2.0f*z);

    f->b[0] = b;
    f->b[1] = b;
    f->a[0] = 1.0f;
    f->a[1] = a;

    return f;
}

// Creates a second order notch filter
// k = Fc/(Fs/2)
IIR_Filter_Res* create_iir_notch_filter(const float k, const float r) {
    const float w0 = PI*k;
    const float a = 2.0f*std::cos(w0);

    auto f = new IIR_Filter_Res(3);
    f->b[0] = 1.0f;
    f->b[1] = -a;
    f->b[2] = 1.0f;

    f->a[0] = 1.0f;
    f->a[1] = -a*r;
    f->a[2] = r*r;

    return f;
}

// Creates a second order peak filter
// k = Fc/(Fs/2)
IIR_Filter_Res* create_iir_peak_filter(const float k, const float r) {
    const float w0 = PI*k;
    const float a = 2.0f*std::cos(w0);

    auto f = new IIR_Filter_Res(3);
    f->b[0] = 0.0f;
    f->b[1] = a*(1.0f-r);
    f->b[2] = r*r - 1.0f;

    f->a[0] = 1.0f;
    f->a[1] = -a*r;
    f->a[2] = r*r;

    return f;
}