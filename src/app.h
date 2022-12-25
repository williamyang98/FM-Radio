#pragma once

#include <stdint.h>
#include <complex>
#include <memory>
#include <thread>

#include "audio/portaudio_output.h"
#include "audio/resampled_pcm_player.h"
#include "audio/portaudio_utility.h"

#include "utility/joint_allocate.h"
#include "utility/observable.h"
#include "utility/reconstruction_buffer.h"

class Broadcast_FM_Demod;
class DifferentialManchesterDecoder;
struct RDS_Decoding_Chain;
struct RDS_Database;

class App 
{
private:
    const int block_size;
    std::unique_ptr<std::thread> runner_thread;

    PaDeviceList pa_devices;
    PortAudio_Output pa_output;
    std::unique_ptr<Resampled_PCM_Player> pcm_player;

    AlignedBlock aligned_block_buf;
    tcb::span<std::complex<uint8_t>> data_u8_buf;
    tcb::span<std::complex<float>> data_f32_buf;
    tcb::span<uint8_t> rds_bytes_decode_buf;

    ReconstructionBuffer<std::complex<uint8_t>> input_buf;
    std::unique_ptr<Broadcast_FM_Demod> broadcast_fm_demod;
    std::unique_ptr<DifferentialManchesterDecoder> differential_manchester_decoder;
    std::unique_ptr<RDS_Decoding_Chain> rds_decoding_chain;

    bool is_output_rds_signal;
    Observable<tcb::span<const float>> obs_on_rds_signal;
public:
    App(const int _block_size);
    ~App();
    size_t Process(tcb::span<const std::complex<uint8_t>> x);
public:
    auto& GetPortAudioOutput() { return pa_output; }
    auto& GetPortAudioDeviceList() { return pa_devices; }
    auto& GetFMDemod() { return *(broadcast_fm_demod.get()); }
    auto& GetIsOutputRDSSignal() { return is_output_rds_signal; }
    RDS_Database& GetRDSDatabase();
private:
    void Run();
};