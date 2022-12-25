#include "rds_group_sync.h"
#include "rds_constants.h"
#include <stdio.h>

#define LOG_MESSAGE(...) fprintf(stderr, "[rds_sync] " __VA_ARGS__)
#define LOG_ERROR(...) fprintf(stderr, "ERROR: [rds_sync] " __VA_ARGS__)

constexpr uint32_t get_bitmask_u32(const int total_bits) {
    return 0xFFFFFFFF >> (32-total_bits);
}

constexpr uint16_t get_bitmask_u16(const int total_bits) {
    return 0xFFFF >> (16-total_bits);
}

// Calculate 10bit CRC
// 26bit -> 16bit data + 10bit crc
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

RDS_Group_Sync::RDS_Group_Sync() {
    rd_block_buf = 0;
    rd_block_buf_bits = 0;

    curr_data_block = 0;
    total_data_block_errors = 0;

    // NOTE: The search for the next block is quite cheap so we don't care
    max_group_desyncs_for_reset = 3;
    curr_groups_desync = 0;
    total_bits_desync = 0;

    state = State::FINDING_SYNC;
}

void RDS_Group_Sync::Process(tcb::span<const uint8_t> x) {
    auto bits = bit_reader_t(x);
    const size_t N = bits.size(); 
    size_t curr_bit = 0;
    while (curr_bit < N) 
    {
        switch (state) {
        case State::FINDING_SYNC:
            curr_bit = FindingSync(bits, curr_bit);
            break;
        case State::READ_BLOCK:
            curr_bit = ReadingGroup(bits, curr_bit);
            break;
        }
    }
}

size_t RDS_Group_Sync::FindingSync(bit_reader_t bits, size_t curr_bit) {
    const size_t N = bits.size();

    // Differential manchester encoding
    // We get the difference between every other symbol
    // It doesn't really matter if we are off by a symbol due to differential encoding
    while (curr_bit < N) {
        PushBit(bits[curr_bit]);
        curr_bit++;

        total_bits_desync++;
        const uint32_t offset_amount = (uint32_t)RDS_BLOCK_OFFSET_VALUES[BlockOffsetID::A];
        const uint32_t input_block = rd_block_buf ^ offset_amount;
        const uint16_t crc_out = CalculateCRC10(input_block);
        if (crc_out != 0) {
            total_bits_desync++;
            continue;
        }

        LOG_MESSAGE("Locked onto block A after %d bits\n", total_bits_desync);
        state = State::READ_BLOCK;
        total_bits_desync = 0;
        rd_block_buf_bits = 0;
        PushBlock(rd_block_buf);
        break;
    }

    return curr_bit;
}

size_t RDS_Group_Sync::ReadingGroup(bit_reader_t bits, size_t curr_bit) {
    const size_t N = bits.size();
    
    // Differential manchester encoding
    // We get the difference between every other symbol
    // It doesn't really matter if we are off by a symbol due to differential encoding
    while (curr_bit < N) {
        PushBit(bits[curr_bit]);
        curr_bit++;

        rd_block_buf_bits++;
        const bool is_block_finished = (rd_block_buf_bits == RDS_PARAMS.nb_block_bits);
        if (!is_block_finished) {
            continue;
        }
        rd_block_buf_bits = 0;

        PushBlock(rd_block_buf);
        if (curr_data_block < RDS_PARAMS.nb_blocks_per_group) {
            continue;
        } 

        obs_on_group.Notify(group);

        const int total_errors = total_data_block_errors;
        curr_data_block = 0;
        total_data_block_errors = 0;
        if (total_errors == 0) {
            curr_groups_desync = 0;
            continue;
        }

        char valid_mask[RDS_PARAMS.nb_blocks_per_group];
        for (int block_id = 0; block_id < RDS_PARAMS.nb_blocks_per_group; block_id++) {
            auto& block = group[block_id];
            valid_mask[block_id] = block.is_valid ? '1' : '0';
        }

        curr_groups_desync++;
        // LOG_MESSAGE("Detected errors in group. valid=%.*s errors=%d (%d/%d)\n", 
        //     RDS_PARAMS.nb_blocks_per_group, valid_mask, total_errors, 
        //     curr_groups_desync, max_group_desyncs_for_reset);

        if (curr_groups_desync >= max_group_desyncs_for_reset) {
            state = State::FINDING_SYNC;
            curr_groups_desync = 0;
            break;
        }
    }

    return curr_bit;
}

void RDS_Group_Sync::PushBit(const bool v) {
    const uint32_t bit = (uint32_t)(v & 0b1);
    constexpr uint32_t BLOCK_BITMASK = get_bitmask_u32(RDS_PARAMS.nb_block_bits);
    rd_block_buf = (rd_block_buf << 1) | bit;
    rd_block_buf = rd_block_buf & BLOCK_BITMASK;
}

// Return number of errors at the end of the data group
void RDS_Group_Sync::PushBlock(const uint32_t x) {
    static auto validate_block = [](const uint32_t x, BlockOffsetID id) -> bool {
        const uint32_t crc_input = x ^ RDS_BLOCK_OFFSET_VALUES[id];
        const uint32_t crc = CalculateCRC10(crc_input);
        return (crc == 0);
    };

    static auto get_data = [](const uint32_t x) -> uint16_t {
        constexpr int nb_data_bits = RDS_PARAMS.nb_block_data_bits;
        constexpr int nb_crc_bits = RDS_PARAMS.nb_block_checksum_bits;
        constexpr uint32_t DATA_MASK = get_bitmask_u32(nb_data_bits) << nb_crc_bits;
        return (x & DATA_MASK) >> nb_crc_bits;
    };

    static auto attempt_decode = [this](const uint32_t x, BlockOffsetID id, rds_block_t& block) -> bool {
        if (!validate_block(x, id)) {
            return false;
        }

        block.block_type = id;
        block.data = get_data(x);
        block.is_valid = true;
        return true;
    };

    // TODO: Add error correction using syndrome

    auto& block = group[curr_data_block];
    block.is_valid = false;

    switch (curr_data_block) {
    case 0: 
        attempt_decode(x, BlockOffsetID::A, block);
        break;
    case 1:
        attempt_decode(x, BlockOffsetID::B, block);
        break;
    case 2:
        attempt_decode(x, BlockOffsetID::C, block) ||
        attempt_decode(x, BlockOffsetID::C1, block);
        break;
    case 3:
        attempt_decode(x, BlockOffsetID::D, block);
        break;
    case 4:
        LOG_ERROR("Invalid group index %d\n", curr_data_block);
        break;
    }
    curr_data_block++;

    if (!block.is_valid) {
        total_data_block_errors++;
    }
}

