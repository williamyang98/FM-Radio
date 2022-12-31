#pragma once

#include <stdint.h>
#include <new>

// Since we need a corresponding aligned deallocator call wrap this in a RAII class
// NOTE: Move and move assign operators are not thread safe 
template <typename T>
class AlignedVector 
{
private:
    T* buf;
    size_t len;
    size_t align;
public:
    AlignedVector(const size_t _len=0, const size_t _align=32) 
    {
        len = _len;
        align = _align;
        buf = NULL;
        if (len > 0) {
            buf = (T*)operator new[](len*sizeof(T), std::align_val_t(align));
        }
    }
    ~AlignedVector() {
        if (buf) {
            operator delete[](buf, std::align_val_t(align));
        }
    }
    AlignedVector(const AlignedVector&) = delete;
    AlignedVector& operator=(const AlignedVector&) = delete;
    // move constructor so we can use RVO
    AlignedVector<T>(AlignedVector<T>&& other) {
        len = other.len;
        align = other.align;
        buf = other.buf;

        other.len = 0;
        other.align = 0;
        other.buf = NULL;
    }
    // move assign operator so we can set by value
    AlignedVector<T>& operator=(AlignedVector<T>&& other) {
        if (buf) {
            operator delete[](buf, std::align_val_t(align));
        }

        len = other.len;
        align = other.align;
        buf = other.buf;

        other.len = 0;
        other.align = 0;
        other.buf = NULL;

        return *this;
    };
    T& operator[](size_t index) {
        return buf[index];
    }
    auto begin() const { return buf; }
    auto end() const { return &buf[len]; }
    auto size() const { return len; }
    auto data() const { return buf; }
};