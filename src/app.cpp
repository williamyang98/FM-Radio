#include "app.h"
#include "demod/broadcast_fm_demod.h"
#include "rds_decoder/rds_init.h"
#include "rds_decoder/differential_manchester_decoder.h"

constexpr size_t SIMD_ALIGN_AMOUNT = 32;

App::App(FILE* _fp_in, FILE* _fp_out, const int _block_size)
: fp_in(_fp_in), fp_out(_fp_out), block_size(_block_size) 
{
    is_running = false;

    aligned_block_buf = AllocateJoint(
        data_u8_buf,            BufferParameters{ (size_t)block_size },
        data_f32_buf,           BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
        rds_bytes_decode_buf,   BufferParameters{ (size_t)16 }
    );

    broadcast_fm_demod = std::make_unique<Broadcast_FM_Demod>(block_size);
    differential_manchester_decoder = std::make_unique<DifferentialManchesterDecoder>(rds_bytes_decode_buf);
    rds_decoding_chain = std::make_unique<RDS_Decoding_Chain>();

    {
        auto& mixer = pa_output.GetMixer();
        auto buf = mixer.CreateManagedBuffer(4);
        auto Fs = pa_output.GetSampleRate();
        pcm_player = std::make_unique<Resampled_PCM_Player>(buf, Fs);

        #ifdef _WIN32
        const auto target_host_api_index = Pa_HostApiTypeIdToHostApiIndex(PORTAUDIO_TARGET_HOST_API_ID);
        const auto target_device_index = Pa_GetHostApiInfo(target_host_api_index)->defaultOutputDevice;
        pa_output.Open(target_device_index);
        #else
        pa_output.Open(Pa_GetDefaultOutputDevice());
        #endif
    }

    broadcast_fm_demod->OnAudioOut().Attach([this](tcb::span<const Frame<float>> x, const int Fs) {
        pcm_player->SetInputSampleRate(Fs);
        pcm_player->ConsumeBuffer(x);
    });
     
    broadcast_fm_demod->OnRDSOut().Attach([this](tcb::span<const float> x) {
        differential_manchester_decoder->Process(x);
        // if (fp_out != NULL) {
        //     fwrite(x.data(), sizeof(float), x.size(), fp_out);
        // }
    });

    differential_manchester_decoder->OnRDSBytes().Attach([this](tcb::span<const uint8_t> x) {
        rds_decoding_chain->Process(x);
    });
}

App::~App() {
    Stop();
}

void App::Start() {
    if (is_running) return;
    is_running = true;
    runner_thread = std::make_unique<std::thread>([this]() { Process(); });
}

void App::Stop() {
    if (!is_running) return;
    is_running = false;
    fclose(fp_in);
    if (runner_thread) {
        runner_thread->join();
    }
}

void App::Process() {
    while (is_running) {
        const size_t nb_read = fread(data_u8_buf.data(), sizeof(std::complex<uint8_t>), block_size, fp_in);
        if (nb_read != block_size) {
            fprintf(stderr, "Failed to read %zu/%zu bytes\n", nb_read, (size_t)block_size);
            break;
        }

        for (int i = 0; i < block_size; i++) {
            data_f32_buf[i] = {
                (float)data_u8_buf[i].real() - 127.0f,
                (float)data_u8_buf[i].imag() - 127.0f,
            };
        }

        broadcast_fm_demod->Process(data_f32_buf);
    }
    is_running = false;
}
    
RDS_Database& App::GetRDSDatabase() {
    return rds_decoding_chain->db;
}
