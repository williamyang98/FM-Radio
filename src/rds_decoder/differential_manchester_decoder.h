#pragma once

#include <stdint.h>
#include "utility/span.h"
#include "utility/observable.h"

class DifferentialManchesterDecoder 
{
private:
    tcb::span<uint8_t> buf;
    size_t curr_byte_index;
    size_t curr_bit_index;
    bool is_read_bit;
    bool prev_bit;

    Observable<tcb::span<const uint8_t>> obs_on_buf;
public: 
    DifferentialManchesterDecoder(tcb::span<uint8_t> _buf) {
        buf = _buf;
        curr_byte_index = 0;
        curr_bit_index = 0;
        is_read_bit = false;
        prev_bit = false;
    }
    void Process(tcb::span<const float> x) {
        for (size_t i = 0; i < x.size(); i++) {
            PushBit(x[i]);
        }
    }
    auto& OnRDSBytes() { return obs_on_buf; }
private:
    void PushBit(float x) {
        // Due to machenster encoding we can skip every other bit
        // Due to the differential nature of the encoding, it doesn't matter which bits
        // we decide to read in, as long as we skip every other bit
        is_read_bit = !is_read_bit;
        if (!is_read_bit) return;

        const bool curr_bit = (x > 0.0f);
        const bool bit = curr_bit ^ prev_bit;
        prev_bit = curr_bit;

        auto& curr_byte = buf[curr_byte_index];
        if (curr_bit_index == 0) {
            curr_byte = 0;
        }

        const uint8_t b = (bit & 0b1);
        curr_byte |= (b << (7-curr_bit_index));

        curr_bit_index++;
        curr_byte_index += (curr_bit_index / 8);
        curr_bit_index   = (curr_bit_index % 8);

        if (curr_byte_index == buf.size()) {
            curr_byte_index = 0;
            obs_on_buf.Notify(buf);
        }

    }
};
