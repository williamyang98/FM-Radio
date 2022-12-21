#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>

#include "broadcast_fm_demod.h"
#include <fftw3.h>
#include "dsp/filter_designer.h"
#include "dsp/calculate_fft.h"
#include "dsp/hilbert_fft_transform.h"
#include "dsp/fftshift.h"
#include "dsp/clamp.h"

#include "demod/bpsk_synchroniser.h"

constexpr size_t SIMD_ALIGN_AMOUNT = 32;

Broadcast_FM_Demod::Broadcast_FM_Demod(const int _block_size)
: block_size(_block_size)
{
    lpf_ds_baseband_factor = 4;     // 256MHz
    lpf_ds_signal_factor = 2;       // 128kHz
    lpf_ds_rds_factor = 8;          // 16kHz
    lpf_ds_audio_stereo_factor = 4; // 32kHz

    Fs_baseband = 1024000;
    Fs_ds_baseband = Fs_baseband/lpf_ds_baseband_factor;    
    Fs_signal = Fs_ds_baseband/lpf_ds_signal_factor;
    Fs_rds = Fs_signal;
    Fs_ds_rds = Fs_rds/lpf_ds_rds_factor;
    Fs_audio_stereo = Fs_signal/lpf_ds_audio_stereo_factor;

    const int block_size_ds_baseband = block_size/lpf_ds_baseband_factor;
    const int block_size_signal = block_size_ds_baseband/lpf_ds_signal_factor;
    const int block_size_ds_rds = block_size_signal/lpf_ds_rds_factor;
    const int block_size_audio_stereo_out = block_size_signal/lpf_ds_audio_stereo_factor;

    fm_demod = std::make_unique<FM_Demod>();

    aligned_block_buf = AllocateJoint(
        // 1. FM demodulation
        lpf_ds_buf,             BufferParameters{ (size_t)block_size_ds_baseband, SIMD_ALIGN_AMOUNT },
        fm_demod_buf,           BufferParameters{ (size_t)block_size_ds_baseband, SIMD_ALIGN_AMOUNT },
        signal_buf,             BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        iq_signal_buf,          BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        // 2. Extract components
        audio_lpr_buf,          BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        audio_lmr_buf,          BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        rds_buf,                BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        pilot_buf,              BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        // 3. PLL
        pll_buf,                BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        pll_lpf_phase_error,    BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        pll_raw_phase_error,    BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        // 4. RDS synchronisation
        rds_ds_buf,             BufferParameters{ (size_t)block_size_ds_rds, SIMD_ALIGN_AMOUNT },
        rds_raw_sym_buf,        BufferParameters{ (size_t)block_size_ds_rds, SIMD_ALIGN_AMOUNT },
        rds_pred_sym_buf,       BufferParameters{ (size_t)block_size_ds_rds, SIMD_ALIGN_AMOUNT },
        // 5. Audio framing
        audio_stereo_mixer_buf, BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        audio_stereo_out_buf,   BufferParameters{ (size_t)block_size_audio_stereo_out, SIMD_ALIGN_AMOUNT },
        // 6. FFT
        fft_baseband_buf,       BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
        fft_mag_baseband_buf,   BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
        fft_lpf_ds_buf,         BufferParameters{ (size_t)block_size_ds_baseband, SIMD_ALIGN_AMOUNT },
        fft_mag_lpf_ds_buf,     BufferParameters{ (size_t)block_size_ds_baseband, SIMD_ALIGN_AMOUNT },
        fft_signal_buf,         BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        fft_mag_signal_buf,     BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        fft_rds_buf,            BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        fft_mag_rds_buf,        BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        fft_audio_lmr_buf,      BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        fft_mag_audio_lmr_buf,  BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        fft_audio_lpr_buf,      BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        fft_mag_audio_lpr_buf,  BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        fft_pilot_buf,          BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        fft_mag_pilot_buf,      BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        fft_pll_pilot_buf,      BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT },
        fft_mag_pll_pilot_buf,  BufferParameters{ (size_t)block_size_signal, SIMD_ALIGN_AMOUNT }
    );
    std::memset(aligned_block_buf.begin(), 0, aligned_block_buf.size());

    // NOTE: We introduce some early rolloff to the downsampling filters
    // This way we make sure we don't alias any higher frequency data into our baseband
    constexpr float DOWNSAMPLING_ROLLOFF_FACTOR = 0.95f;

    // 1. FM demodulation
    // baseband -> [poly_ds_lpf] -> ds_baseband
    if (Fs_baseband != Fs_ds_baseband) {
        const float Fs = (float)Fs_baseband;
        const float Fc = (float)Fs_ds_baseband/2.0f;
        const float k = Fc/(Fs/2.0f) * DOWNSAMPLING_ROLLOFF_FACTOR;
        const int N = filt_cfg.order_poly_ds_lpf_baseband;
        const int M = lpf_ds_baseband_factor;
        const int K = N/M;
        auto* res = create_fir_lpf(k, N-1);
        filt_poly_lpf_ds_baseband = std::make_unique<PolyphaseDownsampler<std::complex<float>>>(res->b, M, K);
        delete res;
    }

    // ds_baseband -> [fm_demodulator] -> [poly_ds_lpf] -> ds_signal
    if (Fs_ds_baseband != Fs_signal) {
        const float Fs = (float)Fs_ds_baseband;
        const float Fc = (float)Fs_signal/2.0f;
        const float k = Fc/(Fs/2.0f) * DOWNSAMPLING_ROLLOFF_FACTOR;
        const int N = filt_cfg.order_poly_ds_lpf_signal;
        const int M = lpf_ds_signal_factor;
        const int K = N/M;
        auto* res = create_fir_lpf(k, N-1);
        filt_poly_lpf_ds_signal = std::make_unique<PolyphaseDownsampler<float>>(res->b, M, K);
        delete res;
    }

    // ds_signal -> [hilbert_transform] -> iq_signal
    {
        const int N = filt_cfg.order_fir_hilbert;
        filt_hilbert_transform = std::make_unique<Hilbert_FIR_Filter>(N);
    }

    // 2. Extract components
    // iq_signal -> [fir_lpf] -> audio_lpr
    {
        const float Fs = (float)Fs_signal;
        const float Fc = (float)params.F_audio_lpr;
        const float k = Fc/(Fs/2.0f);
        // NOTE: Use same order as pilot so PLL is at same delay
        const int N = filt_cfg.order_fir_bpf_pilot;
        auto* res = create_fir_lpf(k, N-1);
        filt_fir_lpf_audio_lpr = std::make_unique<FIR_Filter<std::complex<float>>>(res->b, N);
        delete res;
    }

    // iq_signal -> [fir_bpf] -> audio_lmr
    {
        const float Fs = (float)Fs_signal;
        const float Fcenter = (float)params.F_audio_lmr_center;
        const float Fwidth = (float)params.F_audio_lmr_bandwidth; 
        const float f1 = Fcenter-Fwidth;
        const float f2 = Fcenter+Fwidth;

        const float k1 = f1/(Fs/2.0f);
        const float k2 = f2/(Fs/2.0f);
        // NOTE: Use same order as pilot so PLL is at same delay
        const int N = filt_cfg.order_fir_bpf_pilot;
        auto* res = create_fir_bpf(k1, k2, N-1);
        filt_fir_bpf_audio_lmr = std::make_unique<FIR_Filter<std::complex<float>>>(res->b, N);
        delete res;
    }

    // iq_signal -> [fir_bpf] -> rds
    {
        const float Fs = (float)Fs_signal;
        const float Fcenter = (float)params.F_rds_center;
        const float Fwidth = (float)params.F_rds_bandwidth;
        const float f1 = Fcenter-Fwidth;
        const float f2 = Fcenter+Fwidth;

        const float k1 = f1/(Fs/2.0f);
        const float k2 = f2/(Fs/2.0f);
        // NOTE: Use same order as pilot so PLL is at same delay
        const int N = filt_cfg.order_fir_bpf_pilot;
        auto* res = create_fir_bpf(k1, k2, N-1);
        filt_fir_bpf_rds = std::make_unique<FIR_Filter<std::complex<float>>>(res->b, N);
        delete res;
    }

    // iq_signal -> [fir_bpf] -> pilot
    {
        const float Fs = (float)Fs_signal;
        const float Fc = (float)params.F_pilot;
        const float B = (float)params.F_pilot_deviation;
        const float f1 = Fc-B;
        const float f2 = Fc+B;
        const float k1 = f1/(Fs/2.0f);
        const float k2 = f2/(Fs/2.0f);
        // NOTE: Use same order as pilot so PLL is at same delay
        const int N = filt_cfg.order_fir_bpf_pilot;
        auto* res = create_fir_bpf(k1, k2, N-1);
        filt_fir_bpf_pilot = std::make_unique<FIR_Filter<std::complex<float>>>(res->b, res->N);
        delete res;
    }

    // NOTE: Optional
    {
        const float Fs = (float)Fs_signal;
        const float Fc = (float)params.F_pilot;
        const float k = Fc/(Fs/2.0f);
        const float r = 0.9999f;
        auto* res = create_iir_peak_filter(k, r);
        filt_iir_peak_pilot = std::make_unique<IIR_Filter<std::complex<float>>>(res->b, res->a, res->N);
        delete res;
    }

    // 3. PLL
    // PLL loop filter
    {
        const float Fs = (float)Fs_signal;
        const float Fc = (float)params.F_pilot_deviation;
        const float k = Fc/(Fs/2.0f);
        auto* res = create_iir_single_pole_lpf(k);
        filt_iir_lpf_pll_phase_error = std::make_unique<IIR_Filter<float>>(res->b, res->a, res->N);
        delete res;
    }

    // Setup PLL control loop and fractional mixer
    {
        const float Fs = (float)Fs_signal;
        const float Ts = 1.0f/Fs;
        pll_mixer.f_center = -(float)params.F_pilot; 
        pll_mixer.f_gain = -(float)params.F_pilot_deviation;
        pll_mixer.integrator.KTs = Ts;
        // PI controller 
        pll_prev_phase_error = 0.0f;
        integrator_pll_phase_error.KTs = filt_cfg.pll.integrator_gain*Ts;
    }

    // 4. RDS synchronisation

    // rds_buf -> [poly_ds_lpf] -> ds_rds_buf
    if (Fs_rds != Fs_ds_rds) {
        const float Fs = (float)Fs_rds;
        const float Fc = (float)Fs_ds_rds/2.0f;
        const float k = Fc/(Fs/2.0f) * DOWNSAMPLING_ROLLOFF_FACTOR;
        const int N = filt_cfg.order_poly_ds_lpf_rds;
        const int M = lpf_ds_rds_factor;
        const int K = N/M;
        auto* res = create_fir_lpf(k, N-1);
        filt_poly_lpf_ds_rds = std::make_unique<PolyphaseDownsampler<std::complex<float>>>(res->b, M, K);
        delete res;
    }

    {
        bpsk_sync = std::make_unique<BPSK_Synchroniser>(block_size_ds_rds);

        // AGC so our constellation has normalised points
        // Points: [(0-1j), (0+1j)]
        const auto& cfg = bpsk_sync->GetConfig();
        agc_rds.target_power = cfg.agc_target_power;

        rds_total_symbols = 0;
    }


    // 5. Audio framing
    // NOTE: Optional filter
    // audio_lmr_buf -> [fir_lpf] -> audio_lmr_buf
    {
        const float R = filt_cfg.fir_lpf_lmr_denoise_factor;
        const float Fs = (float)Fs_signal;
        const float Fc = (float)params.F_audio_lmr_bandwidth * R;
        const float k = Fc/(Fs/2.0f);
        const int N = filt_cfg.order_fir_lpf_lmr_denoise;
        auto* res = create_fir_lpf(k, N-1);
        filt_fir_lpf_lmr_denoise = std::make_unique<FIR_Filter<std::complex<float>>>(res->b, N);
        delete res;

        // TODO: Add a delay line to match the L+R component to the L-R filtered component
    }

    // audio_mixer -> [pilot_tone_notch_filter] -> audio_mixer 
    {
        const float Fs = (float)Fs_signal;
        const float Fc = (float)params.F_pilot;
        const float k = Fc/(Fs/2.0f);
        const float r = 0.9999f;
        auto* res = create_iir_notch_filter(k, r);
        filt_iir_notch_pilot = std::make_unique<IIR_Filter<Frame<float>>>(res->b, res->a, res->N);
        delete res;
    }

    // audio_mixer -> [poly_ds_lpf] -> audio_stereo_out
    if (Fs_signal != Fs_audio_stereo) {
        const float Fs = (float)Fs_signal;
        const float Fc = (float)Fs_audio_stereo/2.0f;
        const float k = Fc/(Fs/2.0f) * DOWNSAMPLING_ROLLOFF_FACTOR;
        const int N = filt_cfg.order_poly_ds_lpf_audio;
        const int M = lpf_ds_audio_stereo_factor;
        const int K = N/M;
        auto* res = create_fir_lpf(k, N-1);
        filt_poly_lpf_ds_audio_stereo = std::make_unique<PolyphaseDownsampler<Frame<float>>>(res->b, M, K);
        delete res;
    }

    // 6. FFT
    calc_fft_mag_baseband.SetMode(Calculate_FFT_Mag::Mode::AVERAGE);
    calc_fft_mag_baseband.GetAverageBeta() = 0.1f;
    calc_fft_mag_ds_baseband.SetMode(Calculate_FFT_Mag::Mode::AVERAGE);
    calc_fft_mag_ds_baseband.GetAverageBeta() = 0.1f;
    calc_fft_mag_signal.SetMode(Calculate_FFT_Mag::Mode::AVERAGE);
    calc_fft_mag_signal.GetAverageBeta() = 0.1f;
    calc_fft_mag_rds.SetMode(Calculate_FFT_Mag::Mode::AVERAGE);
    calc_fft_mag_rds.GetAverageBeta() = 0.1f;
    calc_fft_mag_audio_lmr.SetMode(Calculate_FFT_Mag::Mode::AVERAGE);
    calc_fft_mag_audio_lmr.GetAverageBeta() = 0.1f;
    calc_fft_mag_audio_lpr.SetMode(Calculate_FFT_Mag::Mode::AVERAGE);
    calc_fft_mag_audio_lpr.GetAverageBeta() = 0.1f;
    calc_fft_mag_pilot.SetMode(Calculate_FFT_Mag::Mode::AVERAGE);
    calc_fft_mag_pilot.GetAverageBeta() = 0.1f;
    calc_fft_mag_pll_pilot.SetMode(Calculate_FFT_Mag::Mode::AVERAGE);
    calc_fft_mag_pll_pilot.GetAverageBeta() = 0.1f;
}

Broadcast_FM_Demod::~Broadcast_FM_Demod() = default;

void Broadcast_FM_Demod::Process(tcb::span<const std::complex<float>> x)
{
    if (x.size() != block_size) {
        return;
    }

    // Get the FM demodulated signal
    // x -> ... -> signal_buf
    Run_FM_Demodulate(x); 

    // NOTE: Apply Hilbert transform to get IQ version of the FM demodulator output
    //       This way we can perform downmixing using complex multiplication without
    //       having to deal with filtering out image frequencies
    // NOTE: Use the FIR hilbert transform since that doesn't have any edge effects
    //       The FFT version of the hilbert transform works in blocks, which results in
    //       the block boundaries becoming wierd
    //       The FIR hilbert filter is continuous and doesn't suffer from this effect
    // signal -> [hilbert_transform] -> iq_signal
    // HilbertFFTTransform(signal_buf, iq_signal_buf);
    filt_hilbert_transform->process(signal_buf.data(), iq_signal_buf.data(), (int)signal_buf.size());

    // iq_signal_buf -> ... -> [audio_lpr_buf, pilot_buf, audio_lmr_buf, rds_buf]
    FilterSubComponents();

    // pilot_buf -> pll_buf -> [harmonic_mixer] -> [audio_lmr_buf, rds_buf]
    // [audio_lmr_buf, rds_buf] ------^
    ApplyPLL();

    // rds_buf -> ... -> tx_rds_symbols
    SynchroniseRDS();

    // [audio_lpr_buf, audio_lmr_buf] -> ... -> audio_stereo_out_buf
    CombineAudio();

    UpdateFFTs(x);

    // Output demodulated data to listeners
    obs_on_audio_block.Notify(audio_stereo_out_buf, Fs_audio_stereo);
    obs_on_rds_symbols.Notify(rds_pred_sym_buf.first(rds_total_symbols));
}

void Broadcast_FM_Demod::Run_FM_Demodulate(tcb::span<const std::complex<float>> x) {
    // baseband -> [poly_lpf_ds_baseband] -> ds_baseband
    if (filt_poly_lpf_ds_baseband) {
        filt_poly_lpf_ds_baseband->process(x.data(), lpf_ds_buf.data(), (int)lpf_ds_buf.size());
    } else {
        std::copy_n(x.data(), x.size(), lpf_ds_buf.data());
    }

    // ds_baseband -> [fm_demodulator] -> fm_demod
    // TODO: move this to config
    const float Fd = 75e3f;
    const float Fs = (float)Fs_ds_baseband;
    fm_demod->Process(lpf_ds_buf, fm_demod_buf, Fd, Fs);

    // fm_demod -> [poly_lpf_ds_signal] -> signal
    if (filt_poly_lpf_ds_signal) {
        filt_poly_lpf_ds_signal->process(fm_demod_buf.data(), signal_buf.data(), (int)signal_buf.size());
    } else {
        std::copy_n(fm_demod_buf.data(), fm_demod_buf.size(), signal_buf.data());
    }
}

void Broadcast_FM_Demod::FilterSubComponents() {
    // iq_signal -> [fir_audio_lpr] -> audio_lpr
    filt_fir_lpf_audio_lpr->process(iq_signal_buf.data(), audio_lpr_buf.data(), (int)audio_lpr_buf.size());
    // iq_signal -> [fir_bpf] -> audio_lmr
    filt_fir_bpf_audio_lmr->process(iq_signal_buf.data(), audio_lmr_buf.data(), (int)audio_lmr_buf.size());
    // iq_signal -> [fir_bpf] -> rds
    filt_fir_bpf_rds->process(iq_signal_buf.data(), rds_buf.data(), (int)rds_buf.size());

    // Extract pilot tone and normalise to 1 Watt
    if (controls.is_pilot_tone_peak_filter) {
        // iq_signal -> [iir_peak] -> pilot
        // NOTE: Run it through multiple times just to be sure
        filt_iir_peak_pilot->process(iq_signal_buf.data(), pilot_buf.data(), (int)pilot_buf.size());
    } else {
        // iq_signal -> [fir_bpf] -> pilot
        filt_fir_bpf_pilot->process(iq_signal_buf.data(), pilot_buf.data(), (int)pilot_buf.size());
    }
}

void Broadcast_FM_Demod::ApplyPLL() {
    // Guarantee that the pilot signal has the correct power
    // pilot -> [agc] -> pilot
    agc_pilot.process(pilot_buf.data(), pilot_buf.data(), (int)pilot_buf.size());

    // Run our PLL loop
    const size_t N = pilot_buf.size();
    for (size_t i = 0; i < N; i++) {
        // Update phase error with PI controller
        float pll_phase_error_lpf = 0.0f;
        filt_iir_lpf_pll_phase_error->process(&pll_prev_phase_error, &pll_phase_error_lpf, 1);
        integrator_pll_phase_error.process(pll_prev_phase_error);
        integrator_pll_phase_error.yn = clamp(integrator_pll_phase_error.yn, -1.0f, 1.0f);
        const float PI_error = 
            pll_phase_error_lpf*filt_cfg.pll.proportional_gain + 
            integrator_pll_phase_error.yn;
        pll_mixer.phase_error = PI_error;

        // We use the pilot tone as a reference for downmixing the audio L-R stereo and RDS signals
        // This is done by generating fractional harmonics of the PLL pilot tone
        // pll^n = exp(jwt)^n = exp(jnwt);
        const float harmonic_audio_lmr = (float)params.F_audio_lmr_center/(float)params.F_pilot;
        const float harmonic_rds = (float)params.F_rds_center/(float)params.F_pilot;

        const float dt_pilot     = pll_mixer.Update();
        const float dt_audio_lmr = dt_pilot * harmonic_audio_lmr;
        const float dt_rds       = dt_pilot * harmonic_rds;

        static auto get_pll = [](float x) { 
            return std::complex<float>(std::cos(x), std::sin(x)); 
        };

        auto pll_pilot     = get_pll(dt_pilot);
        auto pll_audio_lmr = get_pll(dt_audio_lmr);
        auto pll_rds       = get_pll(dt_rds);

        // Apply fractional PLLs on data components
        audio_lmr_buf[i] *= pll_audio_lmr;
        rds_buf[i] *= pll_rds;

        // Update pilot tone PLL phase error
        // exp(jw0t + phi)*exp(-jw0t) = exp(phi)
        const auto pll_residual = pilot_buf[i] * pll_pilot;
        pll_prev_phase_error = std::atan2(pll_residual.imag(), pll_residual.real());

        pll_buf[i] = pll_pilot;
        pll_raw_phase_error[i] = pll_prev_phase_error;
        pll_lpf_phase_error[i] = PI_error;
    }
}

void Broadcast_FM_Demod::SynchroniseRDS() {
    // rds -> [poly_lpf_ds_rds] -> ds_rds
    if (filt_poly_lpf_ds_rds) {
        filt_poly_lpf_ds_rds->process(rds_buf.data(), rds_ds_buf.data(), (int)rds_ds_buf.size());
    } else {
        std::copy_n(rds_buf.data(), rds_buf.size(), rds_ds_buf.data());
    }

    agc_rds.process(rds_ds_buf.data(), rds_ds_buf.data(), (int)rds_ds_buf.size());
    rds_total_symbols = bpsk_sync->Process(rds_ds_buf, rds_raw_sym_buf);

    for (int i = 0; i < rds_total_symbols; i++) {
        // Output symbols are aligned to the imaginary axis
        const float v = rds_raw_sym_buf[i].imag();
        rds_pred_sym_buf[i] = v;
    }
}

void Broadcast_FM_Demod::CombineAudio() {
    using AudioOut = Broadcast_FM_Demod_Controls::AudioOut;

    // Optional L-R denoiser
    if (controls.is_lmr_lpf && (controls.audio_out != AudioOut::LPR)) {
        filt_fir_lpf_lmr_denoise->process(audio_lmr_buf.data(), audio_lmr_buf.data(), (int)audio_lmr_buf.size());
    }

    // NOTE: The LMR stereo data seems to be on the quadrature component
    // audio_lmr + audio_lpr -> audio_mixer
    const size_t N = audio_stereo_mixer_buf.size();
    switch (controls.audio_out) {
    case AudioOut::STEREO:
        for (size_t i = 0; i < N; i++) {
            const auto lpr = audio_lpr_buf[i].real();
            const auto lmr = audio_lmr_buf[i].imag();
            const float k = controls.audio_stereo_mix_factor;
            audio_stereo_mixer_buf[i] = { 
                lpr+k*lmr, // (L+R) + (L-R) = 2L
                lpr-k*lmr, // (L+R) - (L-R) = 2R
            };
        }
        break;
    case AudioOut::LMR:
        for (size_t i = 0; i < N; i++) {
            const auto lmr = audio_lmr_buf[i].imag();
            audio_stereo_mixer_buf[i] = {lmr, lmr};
        }
        break;
    case AudioOut::LPR:
        for (size_t i = 0; i < N; i++) {
            const auto lpr = audio_lpr_buf[i].real();
            audio_stereo_mixer_buf[i] = {lpr, lpr};
        }
        break;
    }

    // Correct for gain loss due to previous processing
    for (size_t i = 0; i < N; i++) {
        audio_stereo_mixer_buf[i] *= 2.0f;
    }

    // audio_mixer -> [notch_filter] -> audio_mixer
    // Remove the pilot tone which might make it past previous low pass filters
    filt_iir_notch_pilot->process(audio_stereo_mixer_buf.data(), audio_stereo_mixer_buf.data(), (int)audio_stereo_mixer_buf.size());

    // audio_mixer -> [poly_ds_lpf] -> audio_stereo
    if (filt_poly_lpf_ds_audio_stereo) {
        filt_poly_lpf_ds_audio_stereo->process(audio_stereo_mixer_buf.data(), audio_stereo_out_buf.data(), (int)audio_stereo_out_buf.size());
    } else {
        std::copy_n(audio_stereo_mixer_buf.data(), audio_stereo_mixer_buf.size(), audio_stereo_out_buf.data());
    }
}

void Broadcast_FM_Demod::UpdateFFTs(tcb::span<const std::complex<float>> x) {
    CalculateFFT(x, fft_baseband_buf);
    InplaceFFTShift(fft_baseband_buf);
    calc_fft_mag_baseband.Process(fft_baseband_buf, fft_mag_baseband_buf);

    CalculateFFT(lpf_ds_buf, fft_lpf_ds_buf);
    InplaceFFTShift(fft_lpf_ds_buf);
    calc_fft_mag_ds_baseband.Process(fft_lpf_ds_buf, fft_mag_lpf_ds_buf);

    CalculateFFT(iq_signal_buf, fft_signal_buf);
    InplaceFFTShift(fft_signal_buf);
    calc_fft_mag_signal.Process(fft_signal_buf, fft_mag_signal_buf);

    CalculateFFT(rds_buf, fft_rds_buf);
    InplaceFFTShift(fft_rds_buf);
    calc_fft_mag_rds.Process(fft_rds_buf, fft_mag_rds_buf);

    CalculateFFT(audio_lmr_buf, fft_audio_lmr_buf);
    InplaceFFTShift(fft_audio_lmr_buf);
    calc_fft_mag_audio_lmr.Process(fft_audio_lmr_buf, fft_mag_audio_lmr_buf);

    CalculateFFT(audio_lpr_buf, fft_audio_lpr_buf);
    InplaceFFTShift(fft_audio_lpr_buf);
    calc_fft_mag_audio_lpr.Process(fft_audio_lpr_buf, fft_mag_audio_lpr_buf);

    CalculateFFT(pilot_buf, fft_pilot_buf);
    InplaceFFTShift(fft_pilot_buf);
    calc_fft_mag_pilot.Process(fft_pilot_buf, fft_mag_pilot_buf);

    CalculateFFT(pll_buf, fft_pll_pilot_buf);
    InplaceFFTShift(fft_pll_pilot_buf);
    calc_fft_mag_pll_pilot.Process(fft_pll_pilot_buf, fft_mag_pll_pilot_buf);
}