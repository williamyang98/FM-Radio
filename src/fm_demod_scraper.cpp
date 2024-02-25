// Simple SDR app 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <thread>
#include <memory>
#include <complex>
#include <assert.h>
#include <vector>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

#include "audio/frame.h"
#include "fm_demod/broadcast_fm_demod.h"
#include "rds_decoder/differential_manchester_decoder.h"
#include "fm_scraper.h"
#include "getopt/getopt.h"

#include "utility/aligned_allocator.hpp"
#include "utility/joint_allocate.h"
#include "utility/observable.h"
#include "utility/reconstruction_buffer.h"

class App 
{
private:
    const int block_size;

    std::vector<uint8_t, AlignedAllocator<uint8_t>> aligned_block_buf;
    tcb::span<std::complex<uint8_t>> data_u8_buf;
    tcb::span<std::complex<float>> data_f32_buf;
    tcb::span<uint8_t> rds_bytes_decode_buf;

    ReconstructionBuffer<std::complex<uint8_t>> input_buf;
    std::unique_ptr<Broadcast_FM_Demod> broadcast_fm_demod;
    std::unique_ptr<DifferentialManchesterDecoder> differential_manchester_decoder;

    Observable<tcb::span<const uint8_t>> obs_on_rds_bytes;
    Observable<tcb::span<const Frame<float>>, int> obs_on_audio_block;
public:
    explicit App(const int _block_size)
    :   block_size(_block_size),
        input_buf(data_u8_buf) 
    {
        constexpr size_t SIMD_ALIGN_AMOUNT = 32;
        aligned_block_buf = AllocateJoint(
            data_u8_buf,            BufferParameters{ (size_t)block_size },
            data_f32_buf,           BufferParameters{ (size_t)block_size, SIMD_ALIGN_AMOUNT },
            rds_bytes_decode_buf,   BufferParameters{ (size_t)16 }
        );

        broadcast_fm_demod = std::make_unique<Broadcast_FM_Demod>(block_size);
        differential_manchester_decoder = std::make_unique<DifferentialManchesterDecoder>(rds_bytes_decode_buf);

        broadcast_fm_demod->OnAudioOut().Attach([this](tcb::span<const Frame<float>> x, const int Fs) {
            obs_on_audio_block.Notify(x, Fs);
        });
        
        broadcast_fm_demod->OnRDSOut().Attach([this](tcb::span<const float> x) {
            differential_manchester_decoder->Process(x);
        });

        differential_manchester_decoder->OnBytes().Attach([this](tcb::span<const uint8_t> x) {
            obs_on_rds_bytes.Notify(x);
        });
    }
    size_t Process(tcb::span<const std::complex<uint8_t>> x) {
        const size_t N = x.size();
        size_t nb_read = 0;
        while (nb_read < N) {
            nb_read += input_buf.ConsumeBuffer(x.subspan(nb_read));
            if (input_buf.IsFull()) {
                Run();
                input_buf.Reset();
            }
        }
        return nb_read;
    }
public:
    auto& OnAudioBlock() { return obs_on_audio_block; }
    auto& On_RDS_Bytes() { return obs_on_rds_bytes; }
private:
    void Run() {
        for (int i = 0; i < block_size; i++) {
            data_f32_buf[i] = {
                (float)data_u8_buf[i].real() - 127.0f,
                (float)data_u8_buf[i].imag() - 127.0f,
            };
        }
        broadcast_fm_demod->Process(data_f32_buf);
    }
};

struct Arguments {
    const char* input_filename = nullptr;
    const char* output_directory = "data/scraper/default";
    uint32_t block_size = 65536;
};

void usage() {
    const auto args = Arguments();
    assert(args.output_directory != nullptr);

    fprintf(stderr, 
        "fm_demod_scraper, output fm data components to an output directory\n\n"
        "\t[-i input filename (default: %s)]\n"
        "\t    If no file is provided then stdin is used\n"
        "\t[-o output directory (default: %s)]\n"
        "\t[-b block size (default: %u)]\n"
        "\t[-h (show usage)]\n",
        args.input_filename ? args.input_filename : "None",
        args.output_directory,
        args.block_size
    );
}

uint32_t power_ceil(uint32_t x) {
    if (x <= 1) return 1;
    uint32_t power = 2;
    x--;
    while (x >>= 1) power <<= 1;
    return power;
};

Arguments parse_args(int argc, char** argv) {
    Arguments args;
    int block_size = int(args.block_size);

    int opt; 
    while ((opt = getopt_custom(argc, argv, "i:o:b:h")) != -1) {
        switch (opt) {
        case 'i':
            args.input_filename = optarg;
            break;
        case 'o':
            args.output_directory = optarg;
            break;
        case 'b':
            block_size = int(atof(optarg));
            break;
        case 'h':
        default:
            usage();
            exit(0);
        }
    }
    
    if (block_size <= 0) {
        fprintf(stderr, "Block size must be positive (%d)\n", block_size);
        exit(1);
    }

    args.block_size = power_ceil(uint32_t(block_size));
    return args;
}

int main(int argc, char** argv) {
    const auto args = parse_args(argc, argv);

    FILE* fp_in = stdin;
    if (args.input_filename != NULL) {
        fp_in = fopen(args.input_filename, "rb");
        if (fp_in == NULL) {
            fprintf(stderr, "Failed to open file for reading\n");
            return 1;
        }
    }

    assert(args.output_directory != nullptr);

    fprintf(stderr, "Using a block size of %u\n", args.block_size);
#ifdef _WIN32
    _setmode(_fileno(fp_in), _O_BINARY);
#endif

    auto app = App(args.block_size);
    auto fm_scraper = FM_Scraper(args.output_directory);

    app.On_RDS_Bytes().Attach([&fm_scraper](tcb::span<const uint8_t> data) {
        fm_scraper.on_rds_bytes(data);
    });

    app.OnAudioBlock().Attach([&fm_scraper](tcb::span<const Frame<float>> data, const int F_sample) {
        fm_scraper.on_audio_data(data, F_sample);
    });

    // Input loop
    auto input_buf = std::vector<std::complex<uint8_t>>(args.block_size);
    while (true) {
        const size_t nb_read = fread((void*)input_buf.data(), sizeof(std::complex<uint8_t>), input_buf.size(), fp_in);
        if (nb_read != input_buf.size()) {
            fprintf(stderr, "Failed to read %zu/%zu bytes\n", nb_read, input_buf.size());
            break;
        }
        app.Process(input_buf);
    }

    return 0;
}