#define _USE_MATH_DEFINES
#include <cmath>

#include "local_oscillator.h"

LocalOscillator::LocalOscillator() 
{
    dt = 0.0f;
    Ts = 0.0f;
    f_center = 0.0f;
}

std::complex<float> LocalOscillator::Update() {
    const float dx = 2.0f*(float)M_PI*f_center*Ts;
    dt += dx;
    dt = std::fmod(dt, 2.0f*(float)M_PI);
    return std::complex<float>(std::cos(dt), std::sin(dt));
}