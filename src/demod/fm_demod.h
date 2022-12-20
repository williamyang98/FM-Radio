#pragma once

#include <complex>
#include "utility/span.h"

class FM_Demod 
{
private:
    float prev_theta;
public:
    FM_Demod();
    void Reset();
    void Process(
        tcb::span<const std::complex<float>> x, tcb::span<float> y, 
        const float Fd, const float Fs);
};