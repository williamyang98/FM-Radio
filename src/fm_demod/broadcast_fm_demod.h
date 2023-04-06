#pragma once

#include <memory>
#include <complex>

#include "dsp/iir_filter.h"
#include "dsp/fir_filter.h"
#include "dsp/integrator.h"
#include "dsp/agc.h"
#include "dsp/polyphase_filter.h"
#include "dsp/calculate_fft_mag.h"
#include "dsp/hilbert_fir_filter.h"

#include "fm_demod.h"
#include "pll_mixer.h"

#include "audio/frame.h"

#include "utility/joint_allocate.h"
#include "utility/observable.h"

// Forward declare
class BPSK_Synchroniser;

// Fixed parameters of the analogue transmission
struct Broadcast_FM_Demod_Analog_Parameters {
    float F_wbfm_deviation = 75e3f;     // Deviation of WBFM signal
    int F_audio_lpr = 15000;
    int F_pilot = 19000;                // 1st harmonic
    int F_pilot_deviation = 100;
    int F_audio_lmr_center = 38000;     // 2nd harmonic
    int F_audio_lmr_bandwidth = 15000;
    int F_rds_center = 57000;           // 3rd harmonic
    int F_rds_bandwidth = 2000;
    // Number of us of deemphasis filter cutoff 
    // Cutoff of filter starts at fc = 1/(2*pi*T)
    int Tus_min_deemphasis = 1;              
    int Tus_max_deemphasis = 100;              
};

// Fixed demodulator parameters
struct Broadcast_FM_Demod_Config {
    // 1. FM demodulator
    int order_poly_ds_lpf_fm_in = 64;
    int order_poly_ds_lpf_fm_out = 64;
    int order_fir_hilbert = 64+1;         // NOTE: This should be an ODD number for symmetry
    // 2. Lock onto pilot
    struct {
        float integrator_gain = 0.1f;
        float proportional_gain = 0.01f;
    } pilot_pll;
    // 3. Extract components
    int order_poly_ds_lpf_rds = 128;
    int order_poly_ds_lpf_audio = 128;
    float fir_lpf_lmr_denoise_factor = 0.2f; // Amount of bandwidth of L-R to pass with aggressive filtering
    struct {
        float beta_update = 0.1f;            // Correct for phase offset in L-R audio
        int read_stride = 10;
    } audio_lmr_phase;
};

// Controllable demodulator parameters
struct Broadcast_FM_Demod_Controls {
    template <typename T>
    struct EditableControl 
    {
    private:
        bool is_dirty = false;
        T value = 0;
    public:
        void SetValue(T v) { 
            value = v; 
            is_dirty = true;
        }
        void ClearDirty() { is_dirty = false; }
        auto GetValue() const { return value; }
        auto IsDirty() const { return is_dirty; }
    };
    enum AudioOut { LPR, LMR, STEREO };

    AudioOut audio_out = AudioOut::STEREO;
    // Control how much the L-R gets mixed into L+R
    float audio_stereo_mix_factor = 1.0f;
    bool is_use_deemphasis_filter = false;
    EditableControl<int> filt_deemphasis_cutoff;    // us -> fc = 1/(2*pi*T)
    EditableControl<int> filt_audio_lpr_cutoff;     // Hz
    EditableControl<int> filt_audio_lmr_cutoff;     // Hz
};

class Broadcast_FM_Demod 
{
private:
    const int block_size; 
    Broadcast_FM_Demod_Analog_Parameters params;
    Broadcast_FM_Demod_Config filt_cfg;
    Broadcast_FM_Demod_Controls controls;

    // FM contains the following signals
    // Frequency | Description
    // 0-15kHz   | Left+Right mono audio (lpr)
    // 19kHz     | Pilot tone
    // 38±15kHz  | Left-Right stereo difference audio (lmr)
    // 57±2kHz   | RDS digital data (https://en.wikipedia.org/wiki/Radio_Data_System)

    // Demodulator structure
    // 1. FM demodulation
    // baseband -> [poly_ds_lpf] -> fm_in
    // fm_in -> [fm_demodulator] -> fm_demod
    // fm_demod -> [poly_ds_lpf] -> fm_out
    // OPTIONAL: fm_out -> [iir_lpf] -> fm_out
    // fm_out -> [hilbert_transform] -> fm_out_iq

    // 2. Lock onto pilot tone
    // fm_out_iq -> [iir_peak] -> pilot -> [agc] -> pilot

    // 3. Extract components
    // fm_out_iq -> [poly_ds_lpf] -> temp_audio_buf -> [extract_real] -> audio_lpr
    // fm_out_iq -> [pll*2] -> [phase_correct] -> temp_pll_buf -> [poly_ds_lpf] -> temp_audio_buf -> [extract_imag] -> audio_lmr
    // fm_out_iq -> [pll*3] -> temp_pll_buf -> [poly_ds_lpf] -> rds

    // 4. RDS synchronisation
    // rds -> [agc] -> rds
    // rds -> [ted + integrate/dump filter] -> tx_rds_symbols

    // 5. Audio mixing
    // (audio_lpr, audio_lmr) -> [mixer] -> audio_stereo_out

    // All downsampling factors/sample rates in demodulator
    int lpf_ds_fm_in_factor;
    int lpf_ds_fm_out_factor;
    int lpf_ds_rds_factor;
    int lpf_ds_audio_factor;

    int Fs_baseband;
    int Fs_fm_in;
    int Fs_fm_out;
    int Fs_rds;
    int Fs_audio;

    // NOTE: Naming convention of filter
    //       filt_[type of filter]_[output name]
    // 1. FM Demodulation
    std::unique_ptr<PolyphaseDownsampler<std::complex<float>>> filt_poly_ds_lpf_fm_in;
    std::unique_ptr<FM_Demod> fm_demod;
    std::unique_ptr<PolyphaseDownsampler<float>> filt_poly_ds_lpf_fm_out;
    std::unique_ptr<IIR_Filter<float>> filt_iir_lpf_fm_deemphasis;
    std::unique_ptr<Hilbert_FIR_Filter<float>> filt_hilbert_transform;
    // 2. Lock onto pilot tone
    std::unique_ptr<IIR_Filter<std::complex<float>>> filt_iir_peak_pilot;
    std::unique_ptr<IIR_Filter<float>> filt_iir_lpf_pll_phase_error;
    AGC_Filter<std::complex<float>> agc_pilot;
    struct {
        PLL_Mixer mixer;
        Integrator_Block<float> integrator_phase_error;
        float prev_phase_error;
    } pilot_pll;
    // 3. Extract components
    std::unique_ptr<PolyphaseDownsampler<std::complex<float>>> filt_poly_ds_lpf_audio_lpr;
    std::unique_ptr<PolyphaseDownsampler<std::complex<float>>> filt_poly_ds_lpf_audio_lmr;
    std::unique_ptr<PolyphaseDownsampler<std::complex<float>>> filt_poly_ds_lpf_rds;
    float audio_lmr_phase_error;
    // 4. RDS synchronisation
    AGC_Filter<std::complex<float>> agc_rds;
    std::unique_ptr<BPSK_Synchroniser> bpsk_sync;
    int rds_total_symbols;

    // Allocate buffers aligned
    AlignedVector<uint8_t> aligned_block_buf;
    // 1. FM demodulation
    tcb::span<std::complex<float>> fm_in_buf;
    tcb::span<float> fm_demod_buf;
    tcb::span<float> fm_out_buf;
    tcb::span<std::complex<float>> fm_out_iq_buf;
    // 2. Lock onto pilot
    tcb::span<std::complex<float>> pilot_buf;
    tcb::span<float> pll_dt_buf;
    tcb::span<std::complex<float>> pll_buf;
    tcb::span<float> pll_lpf_phase_error;
    tcb::span<float> pll_raw_phase_error;
    // 3. Extract components
    tcb::span<std::complex<float>> temp_pll_buf;
    tcb::span<std::complex<float>> temp_audio_buf;
    tcb::span<float> audio_lmr_buf;
    tcb::span<float> audio_lpr_buf;
    tcb::span<std::complex<float>> rds_buf;
    // 4. RDS synchronisation
    tcb::span<std::complex<float>> rds_raw_sym_buf;
    tcb::span<float> rds_pred_sym_buf;
    // 5. Audio mixing
    tcb::span<Frame<float>> audio_out_buf;
    // 6. FFT
    // 6.1. FM demodulation
    tcb::span<std::complex<float>> fft_baseband_buf;
    tcb::span<float> fft_mag_baseband_buf;
    tcb::span<std::complex<float>> fft_fm_in_buf;
    tcb::span<float> fft_mag_fm_in_buf;
    tcb::span<std::complex<float>> fft_fm_out_buf;
    tcb::span<float> fft_mag_fm_out_buf;
    // 6.2. Lock onto pilot
    tcb::span<std::complex<float>> fft_pilot_buf;
    tcb::span<float> fft_mag_pilot_buf;
    tcb::span<std::complex<float>> fft_pll_pilot_buf;
    tcb::span<float> fft_mag_pll_pilot_buf;
    // 6.3. Extract components
    tcb::span<std::complex<float>> fft_audio_lpr_buf;
    tcb::span<float> fft_mag_audio_lpr_buf;
    tcb::span<std::complex<float>> fft_audio_lmr_buf;
    tcb::span<float> fft_mag_audio_lmr_buf;
    tcb::span<std::complex<float>> fft_rds_buf;
    tcb::span<float> fft_mag_rds_buf;

    // 6.1. FM demodulation
    Calculate_FFT_Mag calc_fft_mag_baseband;
    Calculate_FFT_Mag calc_fft_mag_fm_in;
    Calculate_FFT_Mag calc_fft_mag_fm_out;
    // 6.2. Lock onto pilot
    Calculate_FFT_Mag calc_fft_mag_pilot;
    Calculate_FFT_Mag calc_fft_mag_pll;
    // 6.3. Extract components
    Calculate_FFT_Mag calc_fft_mag_audio_lpr;
    Calculate_FFT_Mag calc_fft_mag_audio_lmr;
    Calculate_FFT_Mag calc_fft_mag_rds;

    // audio_stereo_buf, sample_rate
    Observable<tcb::span<const Frame<float>>, int> obs_on_audio_block;
    Observable<tcb::span<const float>> obs_on_rds_symbols;
public:
    Broadcast_FM_Demod(const int _block_size);
    ~Broadcast_FM_Demod();
    void Process(tcb::span<const std::complex<float>> x);
private:
    void UpdateFilters();
    void Run_FM_Demodulate(tcb::span<const std::complex<float>> x);
    void LockOntoPilot();
    void ExtractComponents();
    void SynchroniseRDS();
    void MixAudio();
public:
    // Internal buffers
    // 1. FM demodulation
    auto GetFMOutIQ() { return fm_out_iq_buf; }
    // 2. Lock onto pilot
    auto GetPilotOutput() { return pilot_buf; }
    auto GetPLLOutput() { return pll_buf; }
    auto Get_PLL_Raw_Phase_Error_Output() { return pll_raw_phase_error; }
    auto Get_PLL_LPF_Phase_Error_Output() { return pll_lpf_phase_error; }
    // 3. Extract components
    auto GetLPRAudioOutput() { return audio_lpr_buf; }
    auto GetLMRAudioOutput() { return audio_lmr_buf; }
    auto GetRDSOutput() { return rds_buf; }
    // 4. RDS synchronisation
    auto GetRDSPredSymbols() { return rds_pred_sym_buf.first(rds_total_symbols); }
    auto GetRDSRawSymbols() { return rds_raw_sym_buf.first(rds_total_symbols); }
    // 5. Audio mixing
    auto GetAudioOut() { return audio_out_buf; }
    // 6. FFT
    // 6.1. FM demodulation
    auto GetBasebandMagnitudeSpectrum() { return fft_mag_baseband_buf; }
    auto GetFMInMagnitudeSpectrum() { return fft_mag_fm_in_buf; }
    auto GetFMOutMagnitudeSpectrum() { return fft_mag_fm_out_buf; }
    // 6.2. Lock onto pilot
    auto GetPilotMagnitudeSpectrum() { return fft_mag_pilot_buf; }
    auto GetPLLPilotMagnitudeSpectrum() { return fft_mag_pll_pilot_buf; }
    // 6.3. Extract components
    auto GetAudioLPRMagnitudeSpectrum() { return fft_mag_audio_lpr_buf; }
    auto GetAudioLMRMagnitudeSpectrum() { return fft_mag_audio_lmr_buf; }
    auto GetRDSMagnitudeSpectrum() { return fft_mag_rds_buf; }

    // Get FFT calculators
    // 6.1. FM demodulation
    auto& GetBasebandMagnitudeSpectrumControls() { return calc_fft_mag_baseband; }
    auto& GetFMInputMagnitudeSpectrumControls() { return calc_fft_mag_fm_in; }
    auto& GetSignalMagnitudeSpectrumControls() { return calc_fft_mag_fm_out; }
    // 6.2. Lock onto pilot
    auto& GetPilotMagnitudeSpectrumControls() { return calc_fft_mag_pilot; }
    auto& GetPLLPilotMagnitudeSpectrumControls() { return calc_fft_mag_pll; }
    // 6.3. Extract components
    auto& GetAudioLPRMagnitudeSpectrumControls() { return calc_fft_mag_audio_lpr; }
    auto& GetAudioLMRMagnitudeSpectrumControls() { return calc_fft_mag_audio_lmr; }
    auto& GetRDSMagnitudeSpectrumControls() { return calc_fft_mag_rds; }

    // Sample rates at various parts of the demodulator
    auto GetBasebandSampleRate() const { return Fs_baseband; }
    auto GetFMInSampleRate() const { return Fs_fm_in; }
    auto GetFMOutSampleRate() const { return Fs_fm_out; }
    auto GetRDSSampleRate() const { return Fs_rds; }
    auto GetAudioSampleRate() const { return Fs_audio; }

    // Internal data/components/settings
    auto GetAudioLMRPhaseError() { return audio_lmr_phase_error; }
    auto& GetBPSKSync() { return *(bpsk_sync.get()); }
    auto& GetAnalogParams() { return params; }
    auto& GetControls() { return controls; }

    // Emit output of demodulator
    auto& OnAudioOut() { return obs_on_audio_block; }
    auto& OnRDSOut() { return obs_on_rds_symbols; }
};