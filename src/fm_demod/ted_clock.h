#pragma once

#include "dsp/integrator.h"

// timing error detector clock
class TED_Clock 
{
public:
    Integrator_Block<float> integrator; // used as ramp oscillator for symbol timing from 0 to 1
    float phase_error;
    float phase_error_gain;      // the phase error is +-1 from zero crossing detector
    float fcenter;
    float fgain;
public:
    TED_Clock();
    // get current voltage in ramp integrator
    float get_current_timing();
    // get a normalised error if oscillator is out of sync
    float get_timing_error();
    // return true if the oscillator resets this sample
    bool update();
};