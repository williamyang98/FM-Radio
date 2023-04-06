#pragma once

#include <memory>
#include <complex>

#include "dsp/integrator.h"
#include "dsp/iir_filter.h"

#include "zero_crossing_detector.h"
#include "trigger_cooldown.h"
#include "ted_clock.h"
#include "pll_mixer.h"

#include "utility/joint_allocate.h"
#include "utility/span.h"

struct BPSK_Synchroniser_Config {
    float F_sample_rate = 16e3f;   
    float F_symbol_rate = 2e3f;
    struct {
        float integrator_gain = 10.0f; 
        float proportional_gain = 0.3f;
    } ted_phase_error;
    struct {
        float integrator_gain = 10.0f;
        float proportional_gain = 0.3f;
    } pll_phase_error;
    float ted_max_freq_offset = 1.5e3f;
    float pll_max_freq_offset = 10.0f;
    float agc_target_power = 0.5f;
};

class BPSK_Synchroniser
{
private:
    const int block_size; 
    BPSK_Synchroniser_Config cfg;

    // Zero crossing detector
    Zero_Crossing_Detector zcd_Q;
    Trigger_Cooldown trigger_cooldown_zcd;

    // TED clock
    TED_Clock ted_clock_int_dump_trig;
    // TED clock PI controller
    Integrator_Block<float> int_ted_phase_error;
    std::unique_ptr<IIR_Filter<float>> filt_iir_lpf_ted_phase_error;
    float ted_prev_phase_error;

    // Integrate and dump trigger
    Integrator_Block<std::complex<float>> int_dump_filter;

    // Symbol PLL to correct for frequency/phase errors
    PLL_Mixer pll_mixer;
    // PLL PI controller
    Integrator_Block<float> int_pll_phase_error;
    std::unique_ptr<IIR_Filter<float>> filt_iir_lpf_pll_phase_error;
    float pll_prev_phase_error;

    AlignedVector<uint8_t> aligned_block_buf;
    // internal buffers
    tcb::span<std::complex<float>> pll_sym_buf;
    tcb::span<bool> zcd_trig_buf;
    tcb::span<bool> int_dump_trigger_buf;
    tcb::span<float> ted_raw_phase_error_buf;
    tcb::span<float> ted_pi_phase_error_buf;
    tcb::span<float> pll_raw_phase_error_buf;
    tcb::span<float> pll_pi_phase_error_buf;
    tcb::span<std::complex<float>> int_dump_filter_buf;
public:
    BPSK_Synchroniser(const int _block_size);
    // NOTE: We expect the input signal to have an average power of 0.5W
    // We output the number of symbols written into y
    int Process(tcb::span<const std::complex<float>> x, tcb::span<std::complex<float>> y);
    const auto& GetConfig() const { return cfg; }
public:
    tcb::span<const std::complex<float>> GetPLLSymbols() const { return pll_sym_buf; }
    tcb::span<const bool> GetZeroCrossings() const { return zcd_trig_buf; }
    tcb::span<const bool> GetIntDumpTriggers() const { return int_dump_trigger_buf; }
    tcb::span<const float> GetTEDRawPhaseError() const { return ted_raw_phase_error_buf; }
    tcb::span<const float> GetTEDPIPhaseError() const { return ted_pi_phase_error_buf; }
    tcb::span<const float> GetPLLRawPhaseError() const { return pll_raw_phase_error_buf; }
    tcb::span<const float> GetPLLPIPhaseError() const { return pll_pi_phase_error_buf; }
    tcb::span<const std::complex<float>> GetIntDumpFilter() const { return int_dump_filter_buf; }
};