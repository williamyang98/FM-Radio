#pragma once

#include <complex>
#include "utility/span.h"

class Calculate_FFT_Mag 
{
public:
    enum Mode { NORMAL, AVERAGE, MAX_HOLD };
private:
    Mode mode;
    float average_beta;
public:
    Calculate_FFT_Mag();
    void Process(tcb::span<const std::complex<float>> x, tcb::span<float> y);
    void SetMode(Mode _mode) { mode = _mode; }
    auto GetMode() const { return mode; }
    float& GetAverageBeta() { return average_beta; }
};