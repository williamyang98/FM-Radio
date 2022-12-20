#include "ted_clock.h"
#include "dsp/clamp.h"

// timing error detector clock
TED_Clock::TED_Clock() {
    phase_error = 0.0f;
    phase_error_gain = 1.0f;      // the phase error is +-1 from zero crossing detector
    fcenter = 0.0f;
    fgain = 0.0f;
}

// get current voltage in ramp integrator
float TED_Clock::get_current_timing() {
    return integrator.yn;
}

// get a normalised error if oscillator is out of sync
float TED_Clock::get_timing_error() {
    float error = get_current_timing();
    error = 2.0f * error;
    // if the zero crossing occurs past the half symbol mark then error is  [-1,0]
    if (error > 1.0f) {
        return error - 2.0f;
    }
    // otherwise if zero crossing occurs before the half symbol mark, then error is [0,1]
    return error;
}

// return true if the oscillator resets this sample
bool TED_Clock::update() {
    float control = phase_error * phase_error_gain;
    control = clamp(control, -1.0f, 1.0f);
    const float freq = fcenter + control*fgain;
    const float v = integrator.process(freq);

    // compensate for finite sampling time
    const float offset = integrator.KTs * freq / 2.0f;
    if (v < (1.0f-offset)) {
        return false;
    }
    // reset integrator otherwise
    integrator.yn = 0.0f;
    return true;
}