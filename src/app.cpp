#include "app.h"
#include "demod/broadcast_fm_demod.h"
#include "rds_decoder/rds_init.h"
#include "rds_decoder/differential_manchester_decoder.h"

constexpr size_t SIMD_ALIGN_AMOUNT = 32;

App::App(const int _block_size)
: block_size(_block_size),
  input_buf(data_u8_buf) 
{
    is_output_rds_signal = false;

    aligned_block_buf = AllocateJoint(
        data_u8_buf,            BufferParameters{ (size_t)block_size },
        data_f32_buf,           BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
        rds_bytes_decode_buf,   BufferParameters{ (size_t)16 }
    );

    broadcast_fm_demod = std::make_unique<Broadcast_FM_Demod>(block_size);
    differential_manchester_decoder = std::make_unique<DifferentialManchesterDecoder>(rds_bytes_decode_buf);
    rds_decoding_chain = std::make_unique<RDS_Decoding_Chain>();

    broadcast_fm_demod->OnAudioOut().Attach([this](tcb::span<const Frame<float>> x, const int Fs) {
        obs_on_audio_block.Notify(x, Fs);
    });
     
    broadcast_fm_demod->OnRDSOut().Attach([this](tcb::span<const float> x) {
        differential_manchester_decoder->Process(x);
        if (is_output_rds_signal) {
            obs_on_rds_signal.Notify(x);
        }
    });

    differential_manchester_decoder->OnRDSBytes().Attach([this](tcb::span<const uint8_t> x) {
        rds_decoding_chain->Process(x);
    });
}

App::~App() = default;

size_t App::Process(tcb::span<const std::complex<uint8_t>> x) {
    const size_t N = x.size();
    size_t nb_read = 0;
    while (nb_read < N) {
        nb_read += input_buf.ConsumeBuffer(x);
        if (input_buf.IsFull()) {
            Run();
            input_buf.Reset();
        }
    }
    return nb_read;
}
    
RDS_Database& App::GetRDSDatabase() {
    return rds_decoding_chain->db;
}

void App::Run() {
    for (int i = 0; i < block_size; i++) {
        data_f32_buf[i] = {
            (float)data_u8_buf[i].real() - 127.0f,
            (float)data_u8_buf[i].imag() - 127.0f,
        };
    }

    broadcast_fm_demod->Process(data_f32_buf);
}
