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

inline static 
std::complex<float> GetPhasor(float dt) { 
    return std::complex<float>(std::cos(dt), std::sin(dt)); 
};

static 
void ConfigureFFTCalc(Calculate_FFT_Mag& calc) {
    calc.SetTrigger(Calculate_FFT_Mag::Trigger::SINGLE);
    calc.SetMode(Calculate_FFT_Mag::Mode::AVERAGE);
    calc.GetAverageBeta() = 0.1f;
};

static 
bool UpdateFFTCalc(
    Calculate_FFT_Mag& calc,
    tcb::span<const std::complex<float>> x, 
    tcb::span<std::complex<float>> y_fft, 
    tcb::span<float> y_fft_mag) 
{
    if (!calc.IsAwaitingUpdate()) {
        return false;
    }
    CalculateFFT(x, y_fft);
    InplaceFFTShift(y_fft);
    calc.Process(y_fft, y_fft_mag);
    return true;
}

template <typename T>
static 
void ProcessPolyrateFilter(PolyphaseDownsampler<T>* filt, tcb::span<const T> x, tcb::span<T> y) 
{
    const size_t N = y.size();
    // We do not have a filter if the input and output sample rates are the same
    if (filt != NULL) {
        filt->process(x.data(), y.data(), (int)N);
    } else {
        assert(x.size() == y.size());
        std::copy_n(x.data(), N, y.data());
    }
}

constexpr 
size_t SIMD_ALIGN_AMOUNT = 32;

Broadcast_FM_Demod::Broadcast_FM_Demod(const int _block_size)
: block_size(_block_size)
{
    lpf_ds_fm_in_factor = 4;    // 256MHz
    lpf_ds_fm_out_factor = 2;   // 128kHz
    lpf_ds_rds_factor = 8;      // 16kHz
    lpf_ds_audio_factor = 4;    // 32kHz

    Fs_baseband = 1024000;
    Fs_fm_in = Fs_baseband/lpf_ds_fm_in_factor;    
    Fs_fm_out = Fs_fm_in/lpf_ds_fm_out_factor;
    Fs_rds = Fs_fm_out/lpf_ds_rds_factor;
    Fs_audio = Fs_fm_out/lpf_ds_audio_factor;

    const int block_size_fm_in = block_size/lpf_ds_fm_in_factor;
    const int block_size_fm_out = block_size_fm_in/lpf_ds_fm_out_factor;
    const int block_size_rds = block_size_fm_out/lpf_ds_rds_factor;
    const int block_size_audio = block_size_fm_out/lpf_ds_audio_factor;

    fm_demod = std::make_unique<FM_Demod>();

    aligned_block_buf = AllocateJoint(
        // 1. FM demodulation
        fm_in_buf,              BufferParameters{ (size_t)block_size_fm_in, SIMD_ALIGN_AMOUNT },
        fm_demod_buf,           BufferParameters{ (size_t)block_size_fm_in, SIMD_ALIGN_AMOUNT },
        fm_out_buf,             BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        fm_out_iq_buf,          BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        // 2. Lock onto pilot
        pilot_buf,              BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        pll_dt_buf,             BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        pll_buf,                BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        pll_lpf_phase_error,    BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        pll_raw_phase_error,    BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        // 3. Extract components
        temp_pll_buf,           BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        temp_audio_buf,         BufferParameters{ (size_t)block_size_audio, SIMD_ALIGN_AMOUNT },
        audio_lpr_buf,          BufferParameters{ (size_t)block_size_audio, SIMD_ALIGN_AMOUNT },
        audio_lmr_buf,          BufferParameters{ (size_t)block_size_audio, SIMD_ALIGN_AMOUNT },
        rds_buf,                BufferParameters{ (size_t)block_size_rds, SIMD_ALIGN_AMOUNT },
        // 4. RDS synchronisation
        rds_raw_sym_buf,        BufferParameters{ (size_t)block_size_rds, SIMD_ALIGN_AMOUNT },
        rds_pred_sym_buf,       BufferParameters{ (size_t)block_size_rds, SIMD_ALIGN_AMOUNT },
        // 5. Audio mixing
        audio_out_buf,          BufferParameters{ (size_t)block_size_audio, SIMD_ALIGN_AMOUNT },
        // 6. FFT
        // 6.1. FM demodulation
        fft_baseband_buf,       BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
        fft_mag_baseband_buf,   BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
        fft_fm_in_buf,          BufferParameters{ (size_t)block_size_fm_in, SIMD_ALIGN_AMOUNT },
        fft_mag_fm_in_buf,      BufferParameters{ (size_t)block_size_fm_in, SIMD_ALIGN_AMOUNT },
        fft_fm_out_buf,         BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        fft_mag_fm_out_buf,     BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        // 6.2. Lock onto pilot
        fft_pilot_buf,          BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        fft_mag_pilot_buf,      BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        fft_pll_pilot_buf,      BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        fft_mag_pll_pilot_buf,  BufferParameters{ (size_t)block_size_fm_out, SIMD_ALIGN_AMOUNT },
        // 6.3. Extract components
        fft_audio_lpr_buf,      BufferParameters{ (size_t)block_size_audio, SIMD_ALIGN_AMOUNT },
        fft_mag_audio_lpr_buf,  BufferParameters{ (size_t)block_size_audio, SIMD_ALIGN_AMOUNT },
        fft_audio_lmr_buf,      BufferParameters{ (size_t)block_size_audio, SIMD_ALIGN_AMOUNT },
        fft_mag_audio_lmr_buf,  BufferParameters{ (size_t)block_size_audio, SIMD_ALIGN_AMOUNT },
        fft_rds_buf,            BufferParameters{ (size_t)block_size_rds, SIMD_ALIGN_AMOUNT },
        fft_mag_rds_buf,        BufferParameters{ (size_t)block_size_rds, SIMD_ALIGN_AMOUNT }
    );
    std::memset(aligned_block_buf.begin(), 0, aligned_block_buf.size());

    // NOTE: We introduce some early rolloff to the downsampling filters
    // This way we make sure we don't alias any higher frequency data into our baseband
    constexpr float DOWNSAMPLING_ROLLOFF_FACTOR = 0.95f;

    // 1. FM demodulation
    // baseband -> [poly_ds_lpf] -> fm_in
    if (Fs_baseband != Fs_fm_in) {
        const float Fs = (float)Fs_baseband;
        const float Fc = (float)Fs_fm_in/2.0f;
        const float k = Fc/(Fs/2.0f) * DOWNSAMPLING_ROLLOFF_FACTOR;
        const int N = filt_cfg.order_poly_ds_lpf_fm_out;
        const int M = lpf_ds_fm_in_factor;
        const int K = N/M;
        auto* res = create_fir_lpf(k, N-1);
        filt_poly_ds_lpf_fm_in = std::make_unique<PolyphaseDownsampler<std::complex<float>>>(res->b, M, K);
        delete res;
    }
    // fm_in -> [fm_demodulator] -> fm_demod -> [poly_ds_lpf] -> fm_out
    if (Fs_fm_in != Fs_fm_out) {
        const float Fs = (float)Fs_fm_in;
        const float Fc = (float)Fs_fm_out/2.0f;
        const float k = Fc/(Fs/2.0f) * DOWNSAMPLING_ROLLOFF_FACTOR;
        const int N = filt_cfg.order_poly_ds_lpf_fm_out;
        const int M = lpf_ds_fm_out_factor;
        const int K = N/M;
        auto* res = create_fir_lpf(k, N-1);
        filt_poly_ds_lpf_fm_out = std::make_unique<PolyphaseDownsampler<float>>(res->b, M, K);
        delete res;
    }
    // fm_out -> [hilbert_transform] -> fm_out_iq
    {
        const int N = filt_cfg.order_fir_hilbert;
        filt_hilbert_transform = std::make_unique<Hilbert_FIR_Filter>(N);
    }


    // 2. Lock onto pilot tone
    // fm_out_iq -> [iir_peak] -> pilot
    {
        const float Fs = (float)Fs_fm_out;
        const float Fc = (float)params.F_pilot;
        const float k = Fc/(Fs/2.0f);
        const float r = 0.9999f;
        auto* res = create_iir_peak_filter(k, r);
        filt_iir_peak_pilot = std::make_unique<IIR_Filter<std::complex<float>>>(res->b, res->a, res->N);
        delete res;
    }
    // PLL loop filter
    {
        const float Fs = (float)Fs_fm_out;
        const float Fc = (float)params.F_pilot_deviation;
        const float k = Fc/(Fs/2.0f);
        auto* res = create_iir_single_pole_lpf(k);
        filt_iir_lpf_pll_phase_error = std::make_unique<IIR_Filter<float>>(res->b, res->a, res->N);
        delete res;
    }
    // Setup PLL control loop and fractional mixer
    {
        const float Fs = (float)Fs_fm_out;
        const float Ts = 1.0f/Fs;
        pilot_pll.mixer.f_center = -(float)params.F_pilot; 
        pilot_pll.mixer.f_gain = -(float)params.F_pilot_deviation;
        pilot_pll.mixer.integrator.KTs = Ts;
        // PI controller 
        pilot_pll.prev_phase_error = 0.0f;
        pilot_pll.integrator_phase_error.KTs = filt_cfg.pilot_pll.integrator_gain*Ts;
    }


    // 3. Extract components
    // fm_out_iq -> [poly_ds_lpf] -> temp_audio -> [extract_real] -> audio_lpr
    if (Fs_fm_out != Fs_audio) {
        const float Fs = (float)Fs_fm_out;
        const float Fc = (float)params.F_audio_lpr;
        const float k = Fc/(Fs/2.0f);
        const int N = filt_cfg.order_poly_ds_lpf_audio;
        const int M = lpf_ds_audio_factor;
        const int K = N/M;
        auto* res = create_fir_lpf(k, N-1);
        filt_poly_ds_lpf_audio_lpr = std::make_unique<PolyphaseDownsampler<std::complex<float>>>(res->b, M, K);
        delete res;
    }
    // fm_out_iq -> [pll*2] -> temp_pll -> [poly_ds_lpf] -> temp_audio -> [phase_correct] -> [extract_real] -> audio_lmr
    if (Fs_fm_out != Fs_audio) {
        // Regular filter
        const float Fs = (float)Fs_fm_out;
        const float Fc = (float)params.F_audio_lmr_bandwidth;
        const float k = Fc/(Fs/2.0f);
        const int N = filt_cfg.order_poly_ds_lpf_audio;
        const int M = lpf_ds_audio_factor;
        const int K = N/M;

        auto* res = create_fir_lpf(k, N-1);
        filt_poly_ds_lpf_audio_lmr = std::make_unique<PolyphaseDownsampler<std::complex<float>>>(res->b, M, K);
        delete res;

        // Aggressive filtering of L-R component if it is noisy
        const float R = filt_cfg.fir_lpf_lmr_denoise_factor;
        const float k_aggressive = k*R;
        auto* res_aggressive = create_fir_lpf(k_aggressive, N-1);
        filt_poly_ds_lpf_audio_lmr_aggressive = std::make_unique<PolyphaseDownsampler<std::complex<float>>>(res_aggressive->b, M, K);
        delete res_aggressive;
    }
    // fm_out_iq -> [pll*3] -> temp_pll -> [poly_ds_lpf] -> rds 
    if (Fs_fm_out != Fs_rds) {
        const float Fs = (float)Fs_fm_out;
        const float Fc = (float)params.F_rds_bandwidth;
        const float k = Fc/(Fs/2.0f);
        const int N = filt_cfg.order_poly_ds_lpf_rds;
        const int M = lpf_ds_rds_factor;
        const int K = N/M;
        auto* res = create_fir_lpf(k, N-1);
        filt_poly_ds_lpf_rds = std::make_unique<PolyphaseDownsampler<std::complex<float>>>(res->b, M, K);
        delete res;
    }
    // phase correction of audio L-R component
    {
        audio_lmr_phase_error = 0.0f;
    }


    // 4. RDS synchronisation
    {
        bpsk_sync = std::make_unique<BPSK_Synchroniser>(block_size_rds);

        // AGC so our constellation has normalised points
        // Points: [(0-1j), (0+1j)]
        const auto& cfg = bpsk_sync->GetConfig();
        agc_rds.target_power = cfg.agc_target_power;

        rds_total_symbols = 0;
    }

    // 6. FFT
    // 6.1. FM demodulation
    ConfigureFFTCalc(calc_fft_mag_baseband);
    ConfigureFFTCalc(calc_fft_mag_fm_in);
    ConfigureFFTCalc(calc_fft_mag_fm_out);
    // 6.2. Lock onto pilot
    ConfigureFFTCalc(calc_fft_mag_pilot);
    ConfigureFFTCalc(calc_fft_mag_pll);
    // 6.3. Extract components
    ConfigureFFTCalc(calc_fft_mag_audio_lpr);
    ConfigureFFTCalc(calc_fft_mag_audio_lmr);
    ConfigureFFTCalc(calc_fft_mag_rds);
}

Broadcast_FM_Demod::~Broadcast_FM_Demod() = default;

void Broadcast_FM_Demod::Process(tcb::span<const std::complex<float>> x)
{
    if (x.size() != block_size) {
        return;
    }

    Run_FM_Demodulate(x); 
    LockOntoPilot();
    ExtractComponents();
    SynchroniseRDS();
    MixAudio();

    // Emit demodulator output
    obs_on_audio_block.Notify(audio_out_buf, Fs_audio);
    obs_on_rds_symbols.Notify(rds_pred_sym_buf.first(rds_total_symbols));
}

void Broadcast_FM_Demod::Run_FM_Demodulate(tcb::span<const std::complex<float>> x) {
    // baseband -> [poly_ds_lpf] -> fm_in
    ProcessPolyrateFilter(filt_poly_ds_lpf_fm_in.get(), x, fm_in_buf);

    // fm_in -> [fm_demodulator] -> fm_demod
    const float Fd = params.F_wbfm_deviation;
    const float Fs = (float)Fs_fm_in;
    fm_demod->Process(fm_in_buf, fm_demod_buf, Fd, Fs);

    // fm_demod -> [poly_ds_lpf] -> fm_out
    ProcessPolyrateFilter<float>(filt_poly_ds_lpf_fm_out.get(), fm_demod_buf, fm_out_buf);

    // fm_out -> [hilbert] -> fm_out_iq
    // HilbertFFTTransform(fm_out_buf, fm_out_iq_buf);
    filt_hilbert_transform->process(fm_out_buf.data(), fm_out_iq_buf.data(), (int)fm_out_buf.size());

    // FFT
    UpdateFFTCalc(calc_fft_mag_baseband, x, fft_baseband_buf, fft_mag_baseband_buf);
    UpdateFFTCalc(calc_fft_mag_fm_in, fm_in_buf, fft_fm_in_buf, fft_mag_fm_in_buf);
    UpdateFFTCalc(calc_fft_mag_fm_out, fm_out_iq_buf, fft_fm_out_buf, fft_mag_fm_out_buf);
}

void Broadcast_FM_Demod::LockOntoPilot() {
    const size_t N = pilot_buf.size();
    // fm_out_iq -> [iir_peak] -> pilot
    filt_iir_peak_pilot->process(fm_out_iq_buf.data(), pilot_buf.data(), (int)N);
    // pilot -> [agc] -> pilot
    agc_pilot.process(pilot_buf.data(), pilot_buf.data(), (int)N);

    // Run our PLL loop
    auto& prev_phase_error = pilot_pll.prev_phase_error;
    auto& integrator_phase_error = pilot_pll.integrator_phase_error;
    auto& mixer = pilot_pll.mixer;
    const auto Kp = filt_cfg.pilot_pll.proportional_gain;
    for (size_t i = 0; i < N; i++) {
        // Update phase error with PI controller
        float phase_error_lpf = 0.0f;
        filt_iir_lpf_pll_phase_error->process(&prev_phase_error, &phase_error_lpf, 1);
        integrator_phase_error.process(prev_phase_error);
        integrator_phase_error.yn = clamp(integrator_phase_error.yn, -1.0f, 1.0f);
        const float PI_error = 
            phase_error_lpf*Kp + 
            integrator_phase_error.yn;
        mixer.phase_error = PI_error;

        const float dt = mixer.Update();
        const auto pll = GetPhasor(dt); 

        // Update pilot tone PLL phase error
        // exp(jw0t + phi)*exp(-jw0t) = exp(phi)
        const auto pll_residual = pilot_buf[i] * pll;
        prev_phase_error = std::atan2(pll_residual.imag(), pll_residual.real());

        pll_dt_buf[i] = dt;
        pll_buf[i] = pll;
        pll_raw_phase_error[i] = prev_phase_error;
        pll_lpf_phase_error[i] = PI_error;
    }

    // FFT
    UpdateFFTCalc(calc_fft_mag_pilot, pilot_buf, fft_pilot_buf, fft_mag_pilot_buf);
    UpdateFFTCalc(calc_fft_mag_pll, pll_buf, fft_pll_pilot_buf, fft_mag_pll_pilot_buf);
}

void Broadcast_FM_Demod::ExtractComponents() {
    // We use the pilot tone as a reference for downmixing the audio L-R stereo and RDS signals
    // This is done by generating fractional harmonics of the PLL pilot tone
    // pll^n = exp(jwt)^n = exp(jnwt);
    const float harmonic_audio_lmr = (float)params.F_audio_lmr_center/(float)params.F_pilot;
    const float harmonic_rds = (float)params.F_rds_center/(float)params.F_pilot;
    const size_t N_pilot = pilot_buf.size();
    const size_t N_audio = audio_lpr_buf.size();
    const size_t N_rds = rds_buf.size();

    // 1. Extract audio L+R component
    // fm_out_iq -> [poly_ds_lpf] -> temp_audio
    ProcessPolyrateFilter<std::complex<float>>(filt_poly_ds_lpf_audio_lpr.get(), fm_out_iq_buf, temp_audio_buf);
    // temp_audio -> [extract_real] -> audio_lpr
    for (size_t i = 0; i < N_audio; i++) {
        audio_lpr_buf[i] = temp_audio_buf[i].real();
    }
    // FFT
    UpdateFFTCalc(calc_fft_mag_audio_lpr, temp_audio_buf, fft_audio_lpr_buf, fft_mag_audio_lpr_buf);

    // 2. Extract audio L-R component
    // fm_out_iq -> [pll*2] -> [phase_correct] -> temp_pll
    for (size_t i = 0; i < N_pilot; i++) {
        // harmonic mixer
        const float dt = pll_dt_buf[i];
        const float dt0 = dt*harmonic_audio_lmr;
        // phase correction
        const float dt1 = dt0 + audio_lmr_phase_error;
        auto pll = GetPhasor(dt1);
        temp_pll_buf[i] = fm_out_iq_buf[i] * pll;
    }
    // temp_pll -> [poly_ds_lpf] -> temp_audio
    if (!controls.is_audio_lmr_aggressive_lpf) {
        ProcessPolyrateFilter<std::complex<float>>(filt_poly_ds_lpf_audio_lmr.get(), temp_pll_buf, temp_audio_buf);
    } else {
        ProcessPolyrateFilter<std::complex<float>>(filt_poly_ds_lpf_audio_lmr_aggressive.get(), temp_pll_buf, temp_audio_buf);
    }
    // Estimate phase error against reference constellation along imaginary axis
    // TODO: We may still end up with a L-R output that is inverted (phase_shift = 180') 
    //       Originally we used a FIR filter so we could estimate the phase delay and 
    //       properly account for the phase shift of the PLL
    //       Now we use an IIR filter for performance reasons
    {
        const size_t stride = (size_t)filt_cfg.audio_lmr_phase.read_stride;
        float avg_phase_error = 0.0f;
        int total_samples = 0;
        for (size_t i = 0; i < N_audio; i+=stride) {
            auto& v = temp_audio_buf[i];
            auto phase = std::atan2(v.imag(), v.real());
            const float estimated_phase_error = 
                (phase > 0.0f) ? 
                (+(float)M_PI/2.0f - phase) :  // +pi/2
                (-(float)M_PI/2.0f - phase);   // -pi/2

            avg_phase_error += estimated_phase_error;
            total_samples++;
        }
        avg_phase_error /= (float)total_samples;

        // Integrate phase error
        const float beta = filt_cfg.audio_lmr_phase.beta_update;
        audio_lmr_phase_error += beta*avg_phase_error;
        audio_lmr_phase_error = std::fmod(audio_lmr_phase_error, 2.0f*(float)M_PI);
    }
    // temp_audio -> [extract_imag] -> audio_lmr
    for (size_t i = 0; i < N_audio; i++) {
        audio_lmr_buf[i] = temp_audio_buf[i].imag();
    }
    // FFT
    UpdateFFTCalc(calc_fft_mag_audio_lmr, temp_audio_buf, fft_audio_lmr_buf, fft_mag_audio_lmr_buf);

    // 3. Extract RDS component
    // fm_out_iq -> [pll*3] -> temp_pll
    for (size_t i = 0; i < N_pilot; i++) {
        const float dt = pll_dt_buf[i];
        const float dt_c = dt * harmonic_rds;
        auto pll = GetPhasor(dt_c);
        temp_pll_buf[i] = fm_out_iq_buf[i] * pll;
    }
    // temp_pll -> [poly_ds_lpf] -> rds
    ProcessPolyrateFilter<std::complex<float>>(filt_poly_ds_lpf_rds.get(), temp_pll_buf, rds_buf);
    // FFT
    UpdateFFTCalc(calc_fft_mag_rds, rds_buf, fft_rds_buf, fft_mag_rds_buf);
}

void Broadcast_FM_Demod::SynchroniseRDS() {
    agc_rds.process(rds_buf.data(), rds_buf.data(), (int)rds_buf.size());
    rds_total_symbols = bpsk_sync->Process(rds_buf, rds_raw_sym_buf);

    for (int i = 0; i < rds_total_symbols; i++) {
        // Output symbols are aligned to the imaginary axis
        const float v = rds_raw_sym_buf[i].imag();
        rds_pred_sym_buf[i] = v;
    }
}

void Broadcast_FM_Demod::MixAudio() {
    using AudioOut = Broadcast_FM_Demod_Controls::AudioOut;

    // NOTE: The LMR stereo data seems to be on the quadrature component
    // (audio_lmr,audio_lpr) -> [mixer] -> audio_out
    const size_t N = audio_out_buf.size();
    switch (controls.audio_out) {
    case AudioOut::STEREO:
        for (size_t i = 0; i < N; i++) {
            const auto lpr = audio_lpr_buf[i];
            const auto lmr = audio_lmr_buf[i];
            const float k = controls.audio_stereo_mix_factor;
            audio_out_buf[i] = { 
                lpr+k*lmr, // (L+R) + (L-R) = 2L
                lpr-k*lmr, // (L+R) - (L-R) = 2R
            };
        }
        break;
    case AudioOut::LMR:
        for (size_t i = 0; i < N; i++) {
            const auto lmr = audio_lmr_buf[i];
            audio_out_buf[i] = {lmr, lmr};
        }
        break;
    case AudioOut::LPR:
        for (size_t i = 0; i < N; i++) {
            const auto lpr = audio_lpr_buf[i];
            audio_out_buf[i] = {lpr, lpr};
        }
        break;
    }

    // Correct for gain loss due to previous processing
    for (size_t i = 0; i < N; i++) {
        audio_out_buf[i] *= 2.0f;
    }
}