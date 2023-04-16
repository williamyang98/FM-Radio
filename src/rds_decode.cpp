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

#include "rds_decoder/rds_decoding_chain.h"
#include "getopt/getopt.h"

struct Arguments {
    const char* input_filename = nullptr;
    uint32_t block_size = 65536;
};

void usage() {
    const auto args = Arguments();

    fprintf(stderr, 
        "rds_decode, decode rds bytes and prints rds data groups to stdout\n\n"
        "\t[-i input filename (default: %s)]\n"
        "\t    If no file is provided then stdin is used\n"
        "\t[-b block size (default: %u)]\n"
        "\t[-h (show usage)]\n",
        args.input_filename ? args.input_filename : "None",
        args.block_size
    );
}

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

    args.block_size = uint32_t(block_size);
    return args;
}

int main(int argc, char** argv) {
    const auto args = parse_args(argc, argv);

    FILE* fp_in = stdin;
    if (args.input_filename != NULL) {
        fp_in = fopen(args.input_filename, "rb");
        if (fp_in == NULL) {
            fprintf(stderr, "Failed to open file for reading (%s)\n", args.input_filename);
            return 1;
        }
    }

    fprintf(stderr, "Using a block size of %u\n", args.block_size);
#ifdef _WIN32
    _setmode(_fileno(fp_in), _O_BINARY);
#endif

    RDS_Decoding_Chain rds_decoding_chain {};

    auto input_buf = std::vector<uint8_t>(args.block_size);
    bool is_reading = true;
    while (is_reading) {
        const size_t nb_read = fread((void*)input_buf.data(), sizeof(uint8_t), input_buf.size(), fp_in);
        if (nb_read != input_buf.size()) {
            fprintf(stderr, "Failed to read %zu/%zu bytes\n", nb_read, input_buf.size());
            is_reading = false;
        }
        rds_decoding_chain.Process(input_buf);
    }

    return 0;
}