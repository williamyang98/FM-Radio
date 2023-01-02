#pragma once

#include "dsp/integrator.h"
#include <complex>

// Creates a local oscillator whose frequency and phase depends on the phase_error
// Fixed parameters: f_center, f_gain, phase_error_gain
// Inputs: phase_error
// Output: Update() -> dt
class PLL_Mixer 
{
public:
    Integrator_Block<float> integrator;
    float phase_error;
    float phase_error_gain;
    float f_center;
    float f_gain;
public:
    PLL_Mixer();
    float Update();
};