#include "calculate_fft_mag.h"
#include <cmath>

Calculate_FFT_Mag::Calculate_FFT_Mag() {
    average_beta = 0.9f;
    mode = Mode::NORMAL;
}

void Calculate_FFT_Mag::Process(
    tcb::span<const std::complex<float>> x, 
    tcb::span<float> y)
{
    const size_t N = y.size();
    const float M = (float)(2*N-1);

    switch (mode) {
    case Mode::NORMAL:
        for (size_t i = 0; i < N; i++) {
            const float v = 20.0f * std::log10(std::abs(x[i]) / M);
            y[i] = v;
        }
        break;
    case Mode::AVERAGE:
        for (size_t i = 0; i < N; i++) {
            const float v = 20.0f * std::log10(std::abs(x[i]) / M);
            const float delta = v-y[i];
            y[i] += average_beta*delta;
        }
        break;
    case Mode::MAX_HOLD:
        for (size_t i = 0; i < N; i++) {
            const float v = 20.0f * std::log10(std::abs(x[i]) / M);
            y[i] = std::max(y[i], v);
        }
        break;
    }
}