#pragma once

#include <vector>
#include <stdio.h>
#include "utility/span.h"

class LoggingBuffer 
{
private:
    std::vector<char> buffer;
    size_t curr_index;
public:
    LoggingBuffer(size_t _max_len=1024) {
        resize(_max_len);
    }
    void resize(size_t len) {
        buffer.resize(len);
        curr_index = 0;
    }

    void reset() { curr_index = 0; }
    size_t capacity() { return buffer.capacity(); }
    size_t length() { return curr_index; }

    template <typename ... U>
    int print(U ... args) {
        char* wr_buffer = &buffer[curr_index];
        const size_t remain_len = capacity()-curr_index;
        const int nb_written = snprintf(wr_buffer, remain_len, args...);
        if (nb_written >= 0) {
            curr_index += (size_t)nb_written;
        }
        return nb_written;
    }

    tcb::span<const char> c_str() {
        return tcb::span(buffer.data(), curr_index);
    }
};