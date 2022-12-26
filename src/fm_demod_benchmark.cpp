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
#include "utility/getopt/getopt.h"

void usage() {
    fprintf(stderr, 
        "fm_demod_benchmark, benchmark the fm demodulator decoding stack\n\n"
        "\t[-i input filename (default: None)]\n"
        "\t    If no file is provided then stdin is used\n"
        "\t[-o output filename (default: None)]\n"
        "\t    If no file is provided then stdout is used\n"
        "\t[-b block size (default: 65536)]\n"
        "\t[-h (show usage)]\n"
    );
}

uint32_t power_ceil(uint32_t x) {
    if (x <= 1) return 1;
    uint32_t power = 2;
    x--;
    while (x >>= 1) power <<= 1;
    return power;
};

int main(int argc, char** argv) {
    uint32_t block_size = 65536;
    const char* rd_filename = NULL;
    const char* wr_filename = NULL;

    int opt; 
    while ((opt = getopt_custom(argc, argv, "i:o:b:h")) != -1) {
        switch (opt) {
        case 'i':
            rd_filename = optarg;
            break;
        case 'o':
            wr_filename = optarg;
            break;
        case 'b':
            block_size = (uint32_t)(atof(optarg));
            break;
        case 'h':
        default:
            usage();
            return 0;
        }
    }

    block_size = power_ceil(block_size);
    if (block_size <= 0) {
        fprintf(stderr, "Block size must be positive (%d)\n", block_size); 
        return 1;
    }

    FILE* fp_in = stdin;
    if (rd_filename != NULL) {
        fp_in = fopen(rd_filename, "rb");
        if (fp_in == NULL) {
            fprintf(stderr, "Failed to open file for reading\n");
            return 1;
        }
    }

    FILE* fp_out = stdout;
    if (wr_filename != NULL) {
        fp_out = fopen(wr_filename, "wb+");
        if (fp_out == NULL) {
            fprintf(stderr, "Failed to open file for writing\n");
            return 1;
        }
    }

    fprintf(stderr, "Using a block size of %u\n", block_size);
#ifdef _WIN32
    _setmode(_fileno(fp_in), _O_BINARY);
    _setmode(_fileno(fp_out), _O_BINARY);
#endif

    // Setup fm demodulator
    auto app = App(block_size);

    // Setup input
    auto input_buf = std::vector<std::complex<uint8_t>>(block_size);
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