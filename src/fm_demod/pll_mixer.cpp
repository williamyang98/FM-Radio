#define _USE_MATH_DEFINES
#include <cmath>

#include "pll_mixer.h"
#include "dsp/clamp.h"

constexpr float PI = (float)M_PI;

PLL_Mixer::PLL_Mixer() {
    phase_error = 0.0f;
    phase_error_gain = 1.0f;
    f_center = 0e3;
    f_gain = 1e3;
}

float PLL_Mixer::Update() {
    float control = phase_error * phase_error_gain;
    control = clamp(control, -1.0f, 1.0f);
    float freq = f_center + control*f_gain;
    float t = integrator.process(2.0f*PI*freq);
    t = std::fmod(t, 2.0f*PI);
    integrator.yn = t;

    return t;
}