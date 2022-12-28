#include "crc10.h"
#include "rds_constants.h"
#include <unordered_map>

constexpr uint16_t get_bitmask_u16(const int total_bits) {
    return 0xFFFF >> (16-total_bits);
}

uint16_t CalculateCRC10(uint32_t x) {
    constexpr int bit_shift = RDS_PARAMS.nb_block_bits-1;
    constexpr uint32_t bit_mask = 1u << bit_shift;
    constexpr uint16_t pop_bit_mask = 1u << RDS_PARAMS.nb_block_checksum_bits;
    constexpr uint16_t crc_bitmask = get_bitmask_u16(RDS_PARAMS.nb_block_checksum_bits);

    uint16_t reg = 0;
    for (int i = 0; i < RDS_PARAMS.nb_block_bits; i++) {
        const uint16_t bit = (uint16_t)((x & bit_mask) >> bit_shift);
        x = x << 1;
        reg = (reg << 1) | bit;
        if (reg & pop_bit_mask) {
            reg = reg ^ RDS_CRC10_POLY;
        }
    }
    return reg & crc_bitmask;
}

// Calculate table for error patterns
auto CRC10_ERROR_PATTERNS = []() {
    auto table = std::unordered_map<uint16_t, uint32_t>();
    const int N = RDS_PARAMS.nb_block_bits;
    const int M = RDS_PARAMS.nb_block_checksum_bits;

    // NOTE: If we are too aggressive with "error correction" we may end
    //       up "recovering" data which passes the CRC but is just meaningless

    // One bit error patterns without editing checksum
    // NOTE: Two bit error patterns gave too many false corrections
    for (int i = M; i < N; i++) {
        const uint32_t error_pattern = (1u << i);
        const uint16_t crc10 = CalculateCRC10(error_pattern);
        table[crc10] = error_pattern;
    }

    // One bit error patterns on checksum
    for (int i = 0; i < M; i++) {
        const uint32_t error_pattern = 1u << i;
        const uint16_t crc10 = CalculateCRC10(error_pattern);
        table[crc10] = error_pattern;
    }

    return std::move(table);
} ();

uint32_t GetCRCErrorFromSyndrome(uint16_t x) {
    auto res = CRC10_ERROR_PATTERNS.find(x);
    if (res == CRC10_ERROR_PATTERNS.end()) {
        return 0;
    }
    return res->second;
}
