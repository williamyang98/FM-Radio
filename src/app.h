#pragma once

#include <stdint.h>
#include <complex>
#include <memory>
#include <thread>

#include "audio/portaudio_output.h"
#include "audio/resampled_pcm_player.h"
#include "audio/portaudio_utility.h"

#include "utility/joint_allocate.h"

class Broadcast_FM_Demod;

class App 
{
private:
    FILE* fp_in;
    FILE* fp_out;
    const int block_size;
    std::unique_ptr<std::thread> runner_thread;
    bool is_running;

    PaDeviceList pa_devices;
    PortAudio_Output pa_output;
    std::unique_ptr<Resampled_PCM_Player> pcm_player;

    AlignedBlock aligned_block_buf;
    tcb::span<std::complex<uint8_t>> data_u8_buf;
    tcb::span<std::complex<float>> data_f32_buf;

    std::unique_ptr<Broadcast_FM_Demod> broadcast_fm_demod;
public:
    App(FILE* _fp_in, FILE* _fp_out, const int _block_size);
    ~App();
public:
    auto& GetPortAudioOutput() { return pa_output; }
    auto& GetPortAudioDeviceList() { return pa_devices; }
    auto& GetFMDemod() { return *(broadcast_fm_demod.get()); }
public:
    void Start();
    void Stop();
private:
    void Process(); 
};