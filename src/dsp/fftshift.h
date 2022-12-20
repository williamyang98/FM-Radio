#pragma once

#include "utility/span.h"

// Shifts FFT results so that the bins match to 0 to Fs
template <typename T>
void FFTShift(tcb::span<const T> x, tcb::span<T> y) {
    const size_t N = x.size();
    const size_t M = N/2;
    // [N/2:N) -> [0:N/2)
    // [0:N/2) -> [N/2:N)
    T tmp;
    for (size_t i = 0, j = M; i < M; i++, j++) {
        y[j] = x[i];
        y[i] = x[j];
    }
}

// Shifts FFT results inplace so that the bins match to 0 to Fs
template <typename T>
void InplaceFFTShift(tcb::span<T> x) {
    const size_t N = x.size();
    const size_t M = N/2;
    // [N/2:N) -> [0:N/2)
    // [0:N/2) -> [N/2:N)
    T tmp;
    for (size_t i = 0, j = M; i < M; i++, j++) {
        tmp = x[i];
        x[i] = x[j];
        x[j] = tmp;
    }
}