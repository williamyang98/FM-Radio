#include "rds_group_sync.h"
#include "rds_constants.h"
#include "crc10.h"

#include <stdio.h>

#define LOG_MESSAGE(...) fprintf(stderr, "[rds_sync] " __VA_ARGS__)
#define LOG_ERROR(...) fprintf(stderr, "ERROR: [rds_sync] " __VA_ARGS__)

constexpr uint32_t get_bitmask_u32(const int total_bits) {
    return 0xFFFFFFFF >> (32-total_bits);
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

// Change data block by reference if we can correct the error
struct ValidateResult {
    bool is_valid = false;
    uint32_t corrected_codeword = 0;
    uint32_t error_pattern = 0;
    uint16_t syndrome = 0;
};

static 
ValidateResult ValidateCRCCodeword(const uint32_t x) {
    auto res = ValidateResult();
    res.corrected_codeword = x;
    res.is_valid = false;
    res.error_pattern = 0;

    res.syndrome = CalculateCRC10(x);
    if (res.syndrome == 0) {
        res.is_valid = true;
        return res;
    }

    res.error_pattern = GetCRCErrorFromSyndrome(res.syndrome);
    // Syndrome has no corresponding error pattern
    if (res.error_pattern == 0) {
        res.is_valid = false;
        return res;
    }

    // We can try correct the error pattern
    const uint32_t x_corr = x ^ res.error_pattern;
    const uint32_t new_syndrome = CalculateCRC10(x_corr);
    if (new_syndrome == 0) {
        res.corrected_codeword = x_corr;
        res.is_valid = true;
        return res;
    }

    res.is_valid = false;
    return res;
};

static 
uint16_t GetDataBits(const uint32_t codeword) {
    constexpr int nb_data_bits = RDS_PARAMS.nb_block_data_bits;
    constexpr int nb_crc_bits = RDS_PARAMS.nb_block_checksum_bits;
    constexpr uint32_t DATA_MASK = get_bitmask_u32(nb_data_bits) << nb_crc_bits;
    return (codeword & DATA_MASK) >> nb_crc_bits;
};

static
const char* GetBlockOffsetName(const BlockOffsetID id) {
    switch (id) {
    case BlockOffsetID::A:  return "A";
    case BlockOffsetID::B:  return "B";
    case BlockOffsetID::C:  return "C";
    case BlockOffsetID::C1: return "C1";
    case BlockOffsetID::D:  return "D";
    case BlockOffsetID::E1: return "E1";
    default:                return "?";
    }
}

static 
bool AttemptDecode(uint32_t x, BlockOffsetID id, rds_block_t& block) {
    x = x ^ RDS_BLOCK_OFFSET_VALUES[id];
    auto res = ValidateCRCCodeword(x);

    if (res.error_pattern != 0) {
        if (res.is_valid) {
            LOG_MESSAGE("Corrected block=%s, error_pattern=%08X\n", 
                GetBlockOffsetName(id), res.error_pattern);
        } else {
            LOG_MESSAGE("Uncorrected block=%s, error_pattern=%08X\n", 
                GetBlockOffsetName(id), res.error_pattern);
        }
    }

    if (!res.is_valid && res.syndrome) {
        LOG_MESSAGE("Uncorrected block=%s, syndrome=%04X\n", 
            GetBlockOffsetName(id), res.syndrome);
    }

    block.block_type = id;
    block.data = GetDataBits(res.corrected_codeword);
    block.is_valid = res.is_valid;
    return block.is_valid;
};

// Return number of errors at the end of the data group
void RDS_Group_Sync::PushBlock(const uint32_t x) {
    auto& block = group[curr_data_block];
    block.is_valid = false;

    switch (curr_data_block) {
    case 0: 
        AttemptDecode(x, BlockOffsetID::A, block);
        break;
    case 1:
        AttemptDecode(x, BlockOffsetID::B, block);
        break;
    case 2:
        AttemptDecode(x, BlockOffsetID::C, block) ||
        AttemptDecode(x, BlockOffsetID::C1, block);
        break;
    case 3:
        AttemptDecode(x, BlockOffsetID::D, block);
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

