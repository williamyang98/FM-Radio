#pragma once

#include <memory>
#include <complex>

#include "dsp/filters.h"
#include "dsp/polyphase_filter.h"
#include "dsp/calculate_fft_mag.h"
#include "dsp/hilbert_fir_filter.h"

#include "demod/fm_demod.h"
#include "demod/pll_mixer.h"

#include "audio/frame.h"

#include "utility/joint_allocate.h"
#include "utility/observable.h"

// Forward declare
class BPSK_Synchroniser;

// Parameters of the analogue transmission
// NOTE: We rely on the pilot signal being the fundemental of the harmonic subcarriers
//       These subcarriers are the stereo data signal and the RDS signal
struct Broadcast_FM_Demod_Analog_Parameters {
    int F_audio_lpr = 15000;
    int F_pilot = 19000;                // 1st harmonic
    int F_pilot_deviation = 100;
    int F_audio_lmr_center = 38000;     // 2nd harmonic
    int F_audio_lmr_bandwidth = 15000;
    int F_rds_center = 57000;           // 3rd harmonic
    int F_rds_bandwidth = 2000;
};

struct Broadcast_FM_Demod_Config {
    int order_poly_ds_lpf_baseband = 32;
    int order_poly_ds_lpf_signal = 32;
    int order_fir_hilbert = 41;         // NOTE: This should be an ODD number for symmetry
    int order_fir_bpf_pilot = 64;       // NOTE: If we are using the BPF for the pilot it 
    int order_fir_lpf_audio_lpr = 64;   //       should have same delay as L+R and L-R
    int order_fir_bpf_audio_lmr = 64; 
    int order_fir_lpf_lmr_denoise = 64;
    int order_poly_ds_lpf_rds = 16;
    int order_poly_ds_lpf_audio = 32;
    float fir_lpf_lmr_denoise_factor = 0.2f; // Amount of bandwidth of L-R to pass
    struct {
        float integrator_gain = 0.1f;
        float proportional_gain = 0.01f;
    } pll;
};

struct Broadcast_FM_Demod_Controls {
    enum AudioOut { LPR, LMR, STEREO };
    AudioOut audio_out = AudioOut::STEREO;
    // NOTE: Control how much the L-R gets mixed into L+R
    //       By default we set this abit low since the stereo data can be abit noisy
    float audio_stereo_mix_factor = 0.5f;
    // NOTE: Higher frequency components in FM are more susceptible to noise
    //       This means the L-R component is quite abit noisier than our L+R component
    //       This is a flag to optionally enable a LPF
    bool is_lmr_lpf = false;
    // NOTE: Change between FIR bpf or IIR peak filter for getting pilot tone
    // NOTE: Using the peak filter means some unknown phase shift in the pilot
    //       This causes the stereo L-R data to be out of phase
    bool is_pilot_tone_peak_filter = false;
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
    // baseband -> [poly_lpf_ds_baseband] -> ds_baseband
    // ds_baseband -> [fm_demodulator] -> fm_demod
    // fm_demod -> [poly_lpf_ds_signal] -> signal
    // signal -> [hilbert_transform] -> iq_signal

    // 2. Extract subcomponents
    // iq_signal -> [lpf_audio_lpr] -> audio_lpr
    // iq_signal -> [bpf_audio_lmr] -> audio_lmr
    // iq_signal -> [bpf_rds] -> rds
    // iq_signal -> [bpf_pilot] -> [agc] -> pilot

    // 3. Apply PLL with fractional mixer
    // pilot -> [pll] -> pll
    // audio_lmr -> [pll] -> audio_lmr
    // rds -> [pll] -> rds

    // 4. RDS synchronisation
    // rds -> [poly_lpf_ds_rds] -> [agc] -> ds_rds
    // ds_rds -> [ted + integrate/dump filter] -> tx_rds_symbols

    // 5. Audio framing
    // audio_lpr + audio_lmr -> [mixer] -> audio_mixer
    // audio_mixer -> [poly_lpf_ds_audio] -> audio_stereo_out

    int lpf_ds_baseband_factor;
    int lpf_ds_signal_factor;
    int lpf_ds_rds_factor;
    int lpf_ds_audio_stereo_factor;

    int Fs_baseband;
    int Fs_ds_baseband;
    int Fs_signal;
    int Fs_rds;
    int Fs_ds_rds;
    int Fs_audio_stereo;

    // 1. FM Demodulation
    std::unique_ptr<PolyphaseDownsampler<std::complex<float>>> filt_poly_lpf_ds_baseband;
    std::unique_ptr<FM_Demod> fm_demod;
    std::unique_ptr<PolyphaseDownsampler<float>> filt_poly_lpf_ds_signal;
    std::unique_ptr<Hilbert_FIR_Filter> filt_hilbert_transform;
    // 2. Extract components
    std::unique_ptr<FIR_Filter<std::complex<float>>> filt_fir_lpf_audio_lpr;
    std::unique_ptr<FIR_Filter<std::complex<float>>> filt_fir_bpf_audio_lmr;
    std::unique_ptr<FIR_Filter<std::complex<float>>> filt_fir_bpf_pilot;
    std::unique_ptr<IIR_Filter<std::complex<float>>> filt_iir_peak_pilot;   // NOTE: We use either the fir_bpf or iir_peak filter
    std::unique_ptr<FIR_Filter<std::complex<float>>> filt_fir_bpf_rds;
    // 3. PLL
    AGC_Filter<std::complex<float>> agc_pilot;
    std::unique_ptr<IIR_Filter<float>> filt_iir_lpf_pll_phase_error;
    PLL_Mixer pll_mixer;
    Integrator_Block<float> integrator_pll_phase_error;
    float pll_prev_phase_error;
    // 4. RDS synchronisation
    std::unique_ptr<PolyphaseDownsampler<std::complex<float>>> filt_poly_lpf_ds_rds;
    AGC_Filter<std::complex<float>> agc_rds;
    std::unique_ptr<BPSK_Synchroniser> bpsk_sync;
    int rds_total_symbols;
    // 5. Audio framing
    // TODO: Change this to a single channel filter
    std::unique_ptr<FIR_Filter<std::complex<float>>> filt_fir_lpf_lmr_denoise;
    std::unique_ptr<IIR_Filter<Frame<float>>> filt_iir_notch_pilot;
    std::unique_ptr<PolyphaseDownsampler<Frame<float>>> filt_poly_lpf_ds_audio_stereo;

    // Allocate buffers aligned
    AlignedBlock aligned_block_buf;
    // 1. FM demodulation
    tcb::span<std::complex<float>> lpf_ds_buf;
    tcb::span<float> fm_demod_buf;
    tcb::span<float> signal_buf;
    tcb::span<std::complex<float>> iq_signal_buf;
    // 2. Extract components
    tcb::span<std::complex<float>> audio_lmr_buf;
    tcb::span<std::complex<float>> audio_lpr_buf;
    tcb::span<std::complex<float>> pilot_buf;
    tcb::span<std::complex<float>> rds_buf;
    // 3. PLL
    tcb::span<std::complex<float>> pll_buf;
    tcb::span<float> pll_lpf_phase_error;
    tcb::span<float> pll_raw_phase_error;
    // 4. RDS synchronisation
    tcb::span<std::complex<float>> rds_ds_buf;
    tcb::span<std::complex<float>> rds_raw_sym_buf;
    tcb::span<float> rds_pred_sym_buf;
    // 5. Audio framing
    tcb::span<Frame<float>> audio_stereo_mixer_buf;
    tcb::span<Frame<float>> audio_stereo_out_buf;
    // 6. FFT
    tcb::span<std::complex<float>> fft_baseband_buf;
    tcb::span<float> fft_mag_baseband_buf;
    tcb::span<std::complex<float>> fft_lpf_ds_buf;
    tcb::span<float> fft_mag_lpf_ds_buf;
    tcb::span<std::complex<float>> fft_signal_buf;
    tcb::span<float> fft_mag_signal_buf;
    tcb::span<std::complex<float>> fft_rds_buf;
    tcb::span<float> fft_mag_rds_buf;
    tcb::span<std::complex<float>> fft_audio_lmr_buf;
    tcb::span<float> fft_mag_audio_lmr_buf;
    tcb::span<std::complex<float>> fft_audio_lpr_buf;
    tcb::span<float> fft_mag_audio_lpr_buf;
    tcb::span<std::complex<float>> fft_pilot_buf;
    tcb::span<float> fft_mag_pilot_buf;
    tcb::span<std::complex<float>> fft_pll_pilot_buf;
    tcb::span<float> fft_mag_pll_pilot_buf;


    Calculate_FFT_Mag calc_fft_mag_baseband;
    Calculate_FFT_Mag calc_fft_mag_ds_baseband;
    Calculate_FFT_Mag calc_fft_mag_signal;
    Calculate_FFT_Mag calc_fft_mag_rds;
    Calculate_FFT_Mag calc_fft_mag_audio_lmr;
    Calculate_FFT_Mag calc_fft_mag_audio_lpr;
    Calculate_FFT_Mag calc_fft_mag_pilot;
    Calculate_FFT_Mag calc_fft_mag_pll_pilot;

    // audio_stereo_buf, sample_rate
    Observable<tcb::span<const Frame<float>>, int> obs_on_audio_block;
    Observable<tcb::span<const float>> obs_on_rds_symbols;
public:
    Broadcast_FM_Demod(const int _block_size);
    ~Broadcast_FM_Demod();
    void Process(tcb::span<const std::complex<float>> x);
private:
    void Run_FM_Demodulate(tcb::span<const std::complex<float>> x);
    void FilterSubComponents();
    void ApplyPLL();
    void SynchroniseRDS();
    void CombineAudio();
    void UpdateFFTs(tcb::span<const std::complex<float>> x);
public:
    // FM demodulation
    tcb::span<const std::complex<float>> GetIQSignal() const { return iq_signal_buf; }
    // Extract components
    tcb::span<const std::complex<float>> GetLPRAudioOutput() const { return audio_lpr_buf; }
    tcb::span<const std::complex<float>> GetLMRAudioOutput() const { return audio_lmr_buf; }
    tcb::span<const std::complex<float>> GetRDSOutput() const { return rds_buf; }
    tcb::span<const std::complex<float>> GetPilotOutput() const { return pilot_buf; }
    // PLL
    tcb::span<const std::complex<float>> GetPLLOutput() const { return pll_buf; }
    tcb::span<const float> Get_PLL_Raw_Phase_Error_Output() const { return pll_raw_phase_error; }
    tcb::span<const float> Get_PLL_LPF_Phase_Error_Output() const { return pll_lpf_phase_error; }
    // RDS synchronisation
    tcb::span<const std::complex<float>> Get_RDS_Downsampled_Output() const { return rds_ds_buf; }
    tcb::span<const float> GetRDSPredSymbols() const { return rds_pred_sym_buf.first(rds_total_symbols); }
    tcb::span<const std::complex<float>> GetRDSRawSymbols() const { return rds_raw_sym_buf.first(rds_total_symbols); }
    auto& GetBPSKSync() { return *(bpsk_sync.get()); }
    // Audio framing
    tcb::span<const Frame<float>> GetStereoOut() const { return audio_stereo_out_buf; }
    // FFT
    tcb::span<const float> GetBasebandMagnitudeSpectrum() const { return fft_mag_baseband_buf; }
    tcb::span<const float> GetDownsampledBasebandMagnitudeSpectrum() const { return fft_mag_lpf_ds_buf; }
    tcb::span<const float> GetSignalMagnitudeSpectrum() const { return fft_mag_signal_buf; }
    tcb::span<const float> GetRDSMagnitudeSpectrum() const { return fft_mag_rds_buf; }
    tcb::span<const float> GetAudioLMRMagnitudeSpectrum() const { return fft_mag_audio_lmr_buf; }
    tcb::span<const float> GetAudioLPRMagnitudeSpectrum() const { return fft_mag_audio_lpr_buf; }
    tcb::span<const float> GetPilotMagnitudeSpectrum() const { return fft_mag_pilot_buf; }
    tcb::span<const float> GetPLLPilotMagnitudeSpectrum() const { return fft_mag_pll_pilot_buf; }
    auto& GetBasebandMagnitudeSpectrumControls() { return calc_fft_mag_baseband; }
    auto& GetDownsampledBasebandMagnitudeSpectrumControls() { return calc_fft_mag_ds_baseband; }
    auto& GetSignalMagnitudeSpectrumControls() { return calc_fft_mag_signal; }
    auto& GetRDSMagnitudeSpectrumControls() { return calc_fft_mag_rds; }
    auto& GetAudioLMRMagnitudeSpectrumControls() { return calc_fft_mag_audio_lmr; }
    auto& GetAudioLPRMagnitudeSpectrumControls() { return calc_fft_mag_audio_lpr; }
    auto& GetPilotMagnitudeSpectrumControls() { return calc_fft_mag_pilot; }
    auto& GetPLLPilotMagnitudeSpectrumControls() { return calc_fft_mag_pll_pilot; }

    auto GetBasebandSampleRate() const { return Fs_baseband; }
    auto GetDownsampledBasebandSampleRate() const { return Fs_ds_baseband; }
    auto GetSignalSampleRate() const { return Fs_signal; }
    auto GetDownsampledRDSSampleRate() const { return Fs_ds_rds; }
    auto GetAudioOutSampleRate() const { return Fs_audio_stereo; }

    auto& GetAnalogParams() { return params; }
    auto& GetControls() { return controls; }

    auto& OnAudioOut() { return obs_on_audio_block; }
    auto& OnRDSOut() { return obs_on_rds_symbols; }
};