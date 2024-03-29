// Simple SDR app 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <thread>
#include <memory>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

#include "app.h"
#include "getopt/getopt.h"

struct Arguments {
    const char* input_filename = nullptr;
    uint32_t block_size = 65536;
};

void usage() {
    const auto args = Arguments();

    fprintf(stderr, 
        "fm_demod_benchmark, benchmark the fm demodulator decoding stack\n\n"
        "\t[-i input filename (default: %s)]\n"
        "\t    If no file is provided then stdin is used\n"
        "\t[-b block size (default: %u)]\n"
        "\t[-h (show usage)]\n",
        args.input_filename ? args.input_filename : "None",
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
    while ((opt = getopt_custom(argc, argv, "i:b:h")) != -1) {
        switch (opt) {
        case 'i':
            args.input_filename = optarg;
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

    fprintf(stderr, "Using a block size of %u\n", args.block_size);
#ifdef _WIN32
    _setmode(_fileno(fp_in), _O_BINARY);
#endif

    // Setup fm demodulator
    auto app = App(args.block_size);

    // Setup input
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