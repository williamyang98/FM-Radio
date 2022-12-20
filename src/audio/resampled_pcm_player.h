#pragma once

#include "frame.h"
#include "ring_buffer.h"
#include <stdint.h>
#include <memory>
#include <vector>

class Resampled_PCM_Player
{
private:
    int input_sample_rate;
    const int output_sample_rate;
    std::shared_ptr<RingBuffer<Frame<float>>> buffer;
    std::vector<Frame<float>> resampling_buffer;
public:
    Resampled_PCM_Player(std::shared_ptr<RingBuffer<Frame<float>>> _buffer, int _output_sample_rate);
    virtual void ConsumeBuffer(tcb::span<const Frame<float>> buf);
    virtual bool SetInputSampleRate(const int _input_sample_rate);
};