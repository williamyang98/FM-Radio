#define _USE_MATH_DEFINES
#include <cmath>

#include "fm_demod.h"

static float wrap_phase(const float x) {
    if (x >= (float)M_PI)       return x - 2.0f*(float)M_PI;
    else if (x <= -(float)M_PI) return x + 2.0f*(float)M_PI;
    else                        return x;
};

FM_Demod::FM_Demod() {
    Reset();
}

void FM_Demod::Reset() {
    prev_theta = 0.0f;
}

// Equation for modulated signal
// m(t) = exp[j*theta(t)]
// theta(t) = integral(w(t), dt)
// w(t) = wc + wd*A(t)

// FM Demodulation
// x1(t) = arg(m(t)) = theta(t) = integral(wc + wd*A(t), dt)
// x2(t) = arg(m(t)*exp(-j*wc*t)) = integral(wd*A(t), dt)
// dx2/dt = (x2[n]-x2[n-1])/Ts = wd*A(t)
// A(t) = (x2[n]-x2[n-1])/(wd*Ts)
void FM_Demod::Process(
    tcb::span<const std::complex<float>> x, tcb::span<float> y, 
    const float Fd, const float Fs)
{
    const size_t N = x.size();
    const float Wd = Fd * 2.0f * (float)M_PI;
    const float Ts = 1.0f / Fs;
    const float A = 1.0f/(Wd*Ts) * 0.5f;

    for (size_t i = 0; i < N; i++) {
        const float curr_theta = std::atan2(x[i].imag(), x[i].real());
        const float delta_theta = wrap_phase(curr_theta-prev_theta);
        y[i] = delta_theta * A;
        prev_theta = curr_theta;
    }
}