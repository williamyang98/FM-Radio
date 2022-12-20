#pragma once

#include <memory>
#include <complex>

#include "dsp/filters.h"
#include "dsp/polyphase_filter.h"
#include "dsp/calculate_fft_mag.h"
#include "dsp/hilbert_fir_filter.h"

#include "demod/fm_demod.h"
#include "demod/pll_mixer.h"
#include "demod/zero_crossing_detector.h"
#include "demod/trigger_cooldown.h"
#include "demod/ted_clock.h"

#include "audio/frame.h"

#include "utility/joint_allocate.h"
#include "utility/observable.h"

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
    int order_fir_bpf_pilot = 64;       // NOTE: This also sets the LMR and RDS bpf orders
    int order_poly_ds_lpf_rds = 16;
    int order_poly_ds_lpf_audio = 32;
    struct {
        float integrator_gain = 1.0f;
        float proportional_gain = 0.1f;
    } pll;
    struct {
        float integrator_gain = 1.0f;
        float proportional_gain = 0.1f;
    } rds_ted;
};

struct Broadcast_FM_Demod_Controls {
    enum AudioOut { LPR, LMR, STEREO };
    AudioOut audio_out = AudioOut::STEREO;
    // NOTE: Control how much the L-R gets mixed into L+R
    //       By default we set this abit low since the stereo data can be abit noisy
    float audio_stereo_mix_factor = 0.5f;
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
    std::unique_ptr<FIR_Filter<std::complex<float>>> filt_fir_bpf_rds;
    AGC_Filter<std::complex<float>> agc_pilot;
    AGC_Filter<std::complex<float>> agc_rds;
    // 3. PLL
    std::unique_ptr<IIR_Filter<float>> filt_iir_lpf_pll_phase_error;
    PLL_Mixer pll_mixer;
    Integrator_Block<float> integrator_pll_phase_error;
    float pll_prev_phase_error;
    // 4. RDS synchronisation
    std::unique_ptr<PolyphaseDownsampler<std::complex<float>>> filt_poly_lpf_ds_rds;
    std::unique_ptr<IIR_Filter<float>> filt_iir_lpf_ted_phase_error;
    Zero_Crossing_Detector zcd_rds;
    Trigger_Cooldown trigger_cooldown_rds_zcd;
    TED_Clock ted_rds_int_dump_trig;
    Integrator_Block<float> integrator_ted_phase_error;
    Integrator_Block<std::complex<float>> rds_int_dump_filter;
    float ted_prev_phase_error;
    Integrator_Block<float> integrator_rds_phase_error;
    // 5. Audio framing
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
    tcb::span<bool> rds_zcd;
    tcb::span<bool> rds_int_dump_trigger;
    tcb::span<float> rds_ted_raw_phase_error;
    tcb::span<float> rds_ted_lpf_phase_error;
    tcb::span<float> rds_pll_phase_error;
    tcb::span<std::complex<float>> rds_int_dump_filter_buf;
    tcb::span<std::complex<float>> rds_bpsk_raw_symbols;
    tcb::span<float> rds_bpsk_symbols;
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

    // NOTE: This is not allocated as part of the joint block buffer
    //       It is instead a subspan of rds_bpsk_symbols which only contains
    //       the symbols demodulated in this frame
    tcb::span<float> tx_rds_symbols;
    tcb::span<std::complex<float>> tx_rds_raw_symbols;

    Calculate_FFT_Mag calc_fft_mag_baseband;
    Calculate_FFT_Mag calc_fft_mag_ds_baseband;
    Calculate_FFT_Mag calc_fft_mag_signal;
    Calculate_FFT_Mag calc_fft_mag_rds;
    Calculate_FFT_Mag calc_fft_mag_audio_lmr;
    Calculate_FFT_Mag calc_fft_mag_audio_lpr;

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
    tcb::span<const bool> Get_RDS_Zero_Crossing() const { return rds_zcd; }
    tcb::span<const bool> Get_RDS_Trig_Dump() const { return rds_int_dump_trigger; }
    tcb::span<const float> Get_RDS_Raw_Phase_Error() const { return rds_ted_raw_phase_error; }
    tcb::span<const float> Get_RDS_LPF_Phase_Error() const { return rds_ted_lpf_phase_error; }
    tcb::span<const std::complex<float>> Get_RDS_Int_Dump_Filter() const { return rds_int_dump_filter_buf; }
    tcb::span<const float> Get_RDS_PLL_Phase_Error() const { return rds_pll_phase_error; }
    tcb::span<const float> Get_TX_RDS_Symbols() const { return tx_rds_symbols; }
    tcb::span<const std::complex<float>> Get_TX_RDS_Raw_Symbols() const { return tx_rds_raw_symbols; }
    // Audio framing
    tcb::span<const Frame<float>> GetStereoOut() const { return audio_stereo_out_buf; }
    // FFT
    tcb::span<const float> GetBasebandMagnitudeSpectrum() const { return fft_mag_baseband_buf; }
    tcb::span<const float> GetDownsampledBasebandMagnitudeSpectrum() const { return fft_mag_lpf_ds_buf; }
    tcb::span<const float> GetSignalMagnitudeSpectrum() const { return fft_mag_signal_buf; }
    tcb::span<const float> GetRDSMagnitudeSpectrum() const { return fft_mag_rds_buf; }
    tcb::span<const float> GetAudioLMRMagnitudeSpectrum() const { return fft_mag_audio_lmr_buf; }
    tcb::span<const float> GetAudioLPRMagnitudeSpectrum() const { return fft_mag_audio_lpr_buf; }
    auto& GetBasebandMagnitudeSpectrumControls() { return calc_fft_mag_baseband; }
    auto& GetDownsampledBasebandMagnitudeSpectrumControls() { return calc_fft_mag_ds_baseband; }
    auto& GetSignalMagnitudeSpectrumControls() { return calc_fft_mag_signal; }
    auto& GetRDSMagnitudeSpectrumControls() { return calc_fft_mag_rds; }
    auto& GetAudioLMRMagnitudeSpectrumControls() { return calc_fft_mag_audio_lmr; }
    auto& GetAudioLPRMagnitudeSpectrumControls() { return calc_fft_mag_audio_lpr; }

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