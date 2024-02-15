#define _USE_MATH_DEFINES
#include <cmath>

#include "bpsk_synchroniser.h"
#include "dsp/filter_designer.h"
#include "dsp/clamp.h"
#include "dsp/simd/chebyshev_sine.h"

constexpr int SIMD_ALIGN_AMOUNT = 32;

BPSK_Synchroniser::BPSK_Synchroniser(const int _block_size) 
: block_size(_block_size)
{
    aligned_block_buf = AllocateJoint(
        pll_sym_buf,            BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
        zcd_trig_buf,           BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
        int_dump_trigger_buf,   BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
        ted_raw_phase_error_buf,BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
        ted_pi_phase_error_buf, BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
        pll_raw_phase_error_buf,BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
        pll_pi_phase_error_buf, BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
        int_dump_filter_buf,    BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT }
    );
    
    // TED loop filter
    {
        const int N = TOTAL_TAPS_IIR_SINGLE_POLE_LPF;
        auto& filt = filt_iir_lpf_ted_phase_error;
        filt = std::make_unique<IIR_Filter<float>>(N);

        const float Fs = (float)cfg.F_sample_rate;
        const float Fc = cfg.ted_max_freq_offset;
        const float k = Fc/(Fs/2.0f);
        create_iir_single_pole_lpf(filt->get_b(), filt->get_a(), k);
    }

    // PLL loop filter
    {
        const int N = TOTAL_TAPS_IIR_SINGLE_POLE_LPF;
        auto& filt = filt_iir_lpf_pll_phase_error;
        filt = std::make_unique<IIR_Filter<float>>(N);

        const float Fs = (float)cfg.F_sample_rate;
        const float Fc = cfg.pll_max_freq_offset;
        const float k = Fc/(Fs/2.0f);
        create_iir_single_pole_lpf(filt->get_b(), filt->get_a(), k);
    }

    // Setup control loop
    {
        const float Fs = cfg.F_sample_rate;
        const float Ts = 1.0f/Fs;
        const float Fsymbol = cfg.F_symbol_rate;
        const int samples_per_symbol = (int)std::round(Fs/Fsymbol);
        const int zcd_cooldown = samples_per_symbol/2;

        // Stop the zero crossing detector from triggering too often
        trigger_cooldown_zcd.N_cooldown = zcd_cooldown;

        // Integrate and dump filter for filtering and acquiring symbols
        // The area under a triangle is given by A = 1/2 * bh
        // b = samples_per_symbol
        // h = 1.0f (normalised bpsk signal)
        const float A = 0.5f * (float)samples_per_symbol * 1.0f;
        int_dump_filter.KTs = 1.0f/A;
        int_dump_filter.yn = 0.0f;

        // Timing clock for triggering the integrate and dump filter
        ted_clock_int_dump_trig.integrator.KTs = Ts;
        ted_clock_int_dump_trig.fcenter = (float)Fsymbol;
        ted_clock_int_dump_trig.fgain = cfg.ted_max_freq_offset;

        // PLL for frequency and phase correction
        pll_mixer.f_center = 0.0f;
        pll_mixer.f_gain = cfg.pll_max_freq_offset;
        pll_mixer.integrator.KTs = Ts;
        pll_mixer.phase_error_gain = 1.0f;

        // NOTE: We scale the update rate of the PI controller based on 
        //       ratio between symbol and sampling rate
        //       This is because the PI controller input is updated every symbol
        const float k = Fsymbol/Fs;
        // PI controller which takes newest phase error and outputs a control signal
        // to the voltage controlled timing clock
        ted_prev_phase_error = 0.0f;
        int_ted_phase_error.KTs = cfg.ted_phase_error.integrator_gain*Ts*k;
        // PI controller for pll mixer
        pll_prev_phase_error = 0.0f;
        int_pll_phase_error.KTs = cfg.pll_phase_error.integrator_gain*Ts*k;
    }
}

int BPSK_Synchroniser::Process(
    tcb::span<const std::complex<float>> x,
    tcb::span<std::complex<float>> y) 
{
    if (x.size() != block_size) {
        return 0;
    }

    // NOTE: We are targeting a BPSK constellation that lies on the imaginary axis
    // Symbol acquisition loop
    int rds_total_symbols = 0;
    for (int i = 0; i < block_size; i++) {
        // PI controller for PLL
        float pll_lpf_phase_error = 0.0f;
        filt_iir_lpf_pll_phase_error->process(&pll_prev_phase_error, &pll_lpf_phase_error, 1);
        int_pll_phase_error.process(pll_prev_phase_error);
        int_pll_phase_error.yn = clamp(int_pll_phase_error.yn, -1.0f, 1.0f);
        const float PI_pll_error = 
            pll_lpf_phase_error*cfg.pll_phase_error.proportional_gain +
            int_pll_phase_error.yn;
        pll_mixer.phase_error = PI_pll_error;

        // Phase correction vector
        // pll = exp(j*theta)
        // v0 = exp(j*phi)
        // v1 = v0 * pll = exp(j*(theta+phi))
        // wrap between [-0.5,+0.5] for chebyshev sine approximation
        const float pll_dt_sin = pll_mixer.Update(); // already wrapped
        float pll_dt_cos = pll_dt_sin+0.25f;
        pll_dt_cos = pll_dt_cos - std::round(pll_dt_cos);
        const auto pll = std::complex<float>(chebyshev_sine(pll_dt_cos), chebyshev_sine(pll_dt_sin)); 
        const auto IQ = x[i] * pll;

        // IQ zero crossing detector on imaginary component
        bool is_zcd = zcd_Q.process(IQ.imag());
        is_zcd = trigger_cooldown_zcd.on_trigger(is_zcd);
        if (is_zcd) {
            ted_prev_phase_error = ted_clock_int_dump_trig.get_timing_error();
        }

        // TED timing clock
        float ted_lpf_phase_error = 0.0f;
        filt_iir_lpf_ted_phase_error->process(&ted_prev_phase_error, &ted_lpf_phase_error, 1);
        // PI controller
        int_ted_phase_error.process(ted_prev_phase_error);
        int_ted_phase_error.yn = clamp(int_ted_phase_error.yn, -1.0f, 1.0f);
        const float PI_ted_error = 
            cfg.ted_phase_error.proportional_gain*ted_lpf_phase_error +
            int_ted_phase_error.yn;
        ted_clock_int_dump_trig.phase_error = -PI_ted_error;

        // Integrate and dump filter
        int_dump_filter.process(IQ);

        // Trigger the integrate and dump filter
        const bool is_ted = ted_clock_int_dump_trig.update();
        if (is_ted) {
            // Dump the integrator output
            const auto sym = int_dump_filter.yn;
            int_dump_filter.yn = {0.0f, 0.0f};


            // Estimate phase error based off known constellation
            // constellation = [(0-1j), (0+1j)]
            // error = [-pi/2, pi/2]
            const float sym_phase = std::atan2(sym.imag(), sym.real());
            const float MAX_PHASE_ERROR = (float)M_PI/2.0f;
            const float estimated_phase_error = 
                (sym_phase > 0.0f) ? 
                (+(float)M_PI/2.0f - sym_phase) :  // +pi/2
                (-(float)M_PI/2.0f - sym_phase);   // -pi/2
            // norm_error = [-1,+1]
            const float normalised_phase_error = estimated_phase_error/MAX_PHASE_ERROR;
            pll_prev_phase_error = normalised_phase_error;

            y[rds_total_symbols] = sym;
            rds_total_symbols++;
        }

        // NOTE: We are keeping track of all the signals for display purposes
        //       An actual implementation would leave all of this out
        pll_sym_buf[i] = IQ;
        zcd_trig_buf[i] = is_zcd;
        int_dump_trigger_buf[i] = is_ted;
        ted_raw_phase_error_buf[i] = ted_prev_phase_error;
        ted_pi_phase_error_buf[i] = PI_ted_error;
        pll_raw_phase_error_buf[i] = pll_prev_phase_error;
        pll_pi_phase_error_buf[i] = PI_pll_error;
        int_dump_filter_buf[i] = int_dump_filter.yn;
    }

    return rds_total_symbols;
}