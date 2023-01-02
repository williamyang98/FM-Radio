#pragma once

#include <cmath>

template <typename T>
class AGC_Filter 
{
public:
    float target_power = 1.0f;
    float current_gain = 0.1f;
    float beta = 0.2f;
    void process(const T* x, T* y, const int N) {
        const float avg_power = calculate_average_power(x, N);
        const float target_gain = std::sqrt(target_power/avg_power);
        current_gain = current_gain + beta*(target_gain - current_gain);
        for (int i = 0; i < N; i++) {
            y[i] = current_gain*x[i];
        }
    }
private:
    float calculate_average_power(const T* x, const int N) {
        float avg_power = 0.0f;
        for (int i = 0; i < N; i++) {
            const float I = x[i].real();
            const float Q = x[i].imag();
            avg_power += (I*I + Q*Q);
        }
        avg_power /= (float)N;
        return avg_power;
    }
};