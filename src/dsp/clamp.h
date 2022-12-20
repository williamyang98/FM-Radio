#pragma once

template <typename T>
T clamp(T x, T _min, T _max) {
    T y = x;
    y = (y > _min) ? y : _min;
    y = (y > _max) ? _max : y;
    return y;
}