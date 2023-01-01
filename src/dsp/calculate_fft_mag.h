#pragma once

#include <complex>
#include "utility/span.h"

class Calculate_FFT_Mag 
{
public:
    enum Mode { NORMAL, AVERAGE, MAX_HOLD };
    enum Trigger { ALWAYS, SINGLE };
private:
    Mode mode;
    Trigger trigger;
    float average_beta;
    bool single_trigger_flag;
public:
    Calculate_FFT_Mag();
    void Process(tcb::span<const std::complex<float>> x, tcb::span<float> y);
    void SetMode(Mode _mode) { mode = _mode; }
    auto GetMode() const { return mode; }
    void SetTrigger(Trigger _trigger) { trigger = _trigger; }
    auto GetTrigger() const { return trigger; }
    void RaiseSingleTrigger() { single_trigger_flag = true; }
    float& GetAverageBeta() { return average_beta; }
    bool IsAwaitingUpdate() {
        if (trigger == Trigger::ALWAYS) return true;
        return single_trigger_flag;
    }
};