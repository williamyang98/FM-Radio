#pragma once

#include <stdint.h>
#include "utility/span.h"
#include "utility/observable.h"
#include "rds_constants.h"

class RDS_Group_Sync 
{
private:
    // Read bits into block buffer
    uint32_t rd_block_buf;
    int rd_block_buf_bits;

    class bit_reader_t 
    {
    private:         
        tcb::span<const uint8_t> x;
        const size_t nb_bits;
    public:
        bit_reader_t(tcb::span<const uint8_t> _x)
        : x(_x), nb_bits(_x.size()*8) {}
        bool operator[](const size_t curr_bit) {
            const size_t curr_byte = curr_bit / 8; 
            const size_t bit_idx   = curr_bit % 8;
            const size_t shift = 7-bit_idx;
            return (x[curr_byte] & (1 << shift)) >> shift;
        }
        size_t size() { return nb_bits; }
    };

    // Process 26bit data+crc into 16bit data
    rds_group_t group;
    int curr_data_block;
    int total_data_block_errors;

    // Check if we should resync to block A
    int max_group_desyncs_for_reset;
    int curr_groups_desync;
    int total_bits_desync;

    enum State { FINDING_SYNC, READ_BLOCK };
    State state;

    Observable<rds_group_t> obs_on_group;
public:
    RDS_Group_Sync();
    void Process(tcb::span<const uint8_t> x);
    auto& OnGroup() { return obs_on_group; }
private:
    size_t FindingSync(bit_reader_t bits, size_t curr_bit);
    size_t ReadingGroup(bit_reader_t bits, size_t curr_bit);
    void PushBit(const bool v);
    void PushBlock(const uint32_t x);
};