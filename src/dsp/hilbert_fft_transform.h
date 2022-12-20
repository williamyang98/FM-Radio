#pragma once

#include "calculate_fft.h"
#include "utility/span.h"

// Performs the Hilbert transform using spectral manipulation via FFT
void HilbertFFTTransform(tcb::span<const float> x, tcb::span<std::complex<float>> y) {
    const size_t N = x.size();
    const size_t M = N/2;
    for (size_t i = 0; i < N; i++) {
        y[i] = { x[i], 0.0f };
    }

    CalculateFFT(y, y);

    // Set negative frequencies to 0
    for (size_t i = M; i < N; i++) {
        y[i] = { 0.0f, 0.0f };
    }

    CalculateIFFT(y, y);

    // Rescale due to removal of half the spectrum
    const float A = 1.0f/(float)N;
    for (size_t i = 0; i < N; i++) {
        y[i] *= A;
    }
}
