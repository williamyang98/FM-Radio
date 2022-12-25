#pragma once
#include <stdint.h>

// Clause 2.1: Baseband coding structure
constexpr struct {
    int nb_block_bits = 26;
    int nb_block_data_bits = 16;
    int nb_block_checksum_bits = 10;
    int nb_blocks_per_group = 4;
} RDS_PARAMS;

// Clause 2.3: Error protection
// g(x) = x^10 + x^8 + x^7 + x^5 + x^4 + x^3 + 1
constexpr uint16_t RDS_CRC10_POLY = 0b0110111001;
// constexpr uint16_t RDS_CRC10_POLY = 0b0011101101; // reversed

// Annex A: Offset words to be used for group and block synchronisation
// Table A.1: Offset words d(x)
constexpr int TOTAL_BLOCK_OFFSET_VALUES = 6;
constexpr uint16_t RDS_BLOCK_OFFSET_VALUES[TOTAL_BLOCK_OFFSET_VALUES] = {
    0b0011111100, // A  
    0b0110011000, // B  
    0b0101101000, // C  
    0b1101010000, // C1   
    0b0110110100, // D  
    0b0000000000, // E1  
};
enum BlockOffsetID: int { A=0, B=1, C=2, C1=3, D=4, E1=5 };

struct rds_block_t {
    uint16_t data;
    BlockOffsetID block_type;
    bool is_valid;
};

struct rds_group_t {
    rds_block_t blocks[RDS_PARAMS.nb_blocks_per_group];
    auto& operator[](size_t i) { return blocks[i]; }
};